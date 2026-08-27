#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
export PROJECT

source "${SCRIPT_DIR}/tints_config.env"

MODE="${1:---plan}"
case "$MODE" in
    --plan|--run) ;;
    *)
        echo "Usage: $0 [--plan|--run]" >&2
        exit 2
        ;;
esac

ANA="${PROJECT}/06_genomics"
QC_TABLE="${ANA}/07_checkm/library_taxonomy_qc.tsv"

ROOT="${ANA}/07R_salmonella_rescue"
CANDIDATES="${ROOT}/rescue_candidates.tsv"
SUMMARY="${ROOT}/rescue_summary.tsv"
PROV="${ROOT}/rescue_provenance.tsv"

# Rescue TRIAGE criteria.
# These do NOT replace final genome QC.
RESCUE_MIN_SALMONELLA_PCT="20"
RESCUE_MIN_DEPTH_X="30"
RESCUE_CAP_X="120"
RESCUE_GENOME_SIZE_BP="4800000"

KRAKENTOOLS="${PROJECT}/tools/KrakenTools/extract_kraken_reads.py"
RESCUE_PY="${HOME}/miniforge3/envs/tints-rescue/bin/python"

mkdir -p "$ROOT"

for x in \
    "$QC_TABLE" \
    "$KRAKENTOOLS" \
    "$RESCUE_PY" \
    "$KRAKEN2_BIN" \
    "$RASUSA_BIN" \
    "$SPADES_BIN" \
    "$CONDA" \
    "$MLST_BIN" \
    "$ANY2FASTA_BIN" \
    "$KRAKEN_DB" \
    "$AMRFINDER_DB" \
    "$NCBI_BLAST_BIN_DIR"
do
    [[ -e "$x" ]] || {
        echo "ERROR: required path missing: $x" >&2
        exit 1
    }
done

{
    printf 'key\tvalue\n'
    printf 'created_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'mode\t%s\n' "$MODE"
    printf 'runner_sha256\t%s\n' "$(shasum -a 256 "$0" | awk '{print $1}')"
    printf 'rescue_min_salmonella_pct\t%s\n' "$RESCUE_MIN_SALMONELLA_PCT"
    printf 'rescue_min_depth_x\t%s\n' "$RESCUE_MIN_DEPTH_X"
    printf 'rescue_cap_x\t%s\n' "$RESCUE_CAP_X"
    printf 'rescue_genome_size_bp\t%s\n' "$RESCUE_GENOME_SIZE_BP"
    printf 'kraken_db\t%s\n' "$KRAKEN_DB"
    printf 'kraken_db_snapshot\t%s\n' "$KRAKEN_DB_SNAPSHOT"
    printf 'kraken_manifest_sha256\t%s\n' "$KRAKEN_DB_MANIFEST_SHA256"
    printf 'krakentools_commit\t%s\n' \
        "$(git -C "${PROJECT}/tools/KrakenTools" rev-parse HEAD)"
    printf 'spades_version\t%s\n' \
        "$("$SPADES_BIN" --version 2>&1 | head -1)"
    printf 'rasusa_version\t%s\n' \
        "$("$RASUSA_BIN" --version 2>&1 | head -1)"
    printf 'checkm_env\t%s\n' "$CHECKM_ENV"
    printf 'quast_env\t%s\n' "$QUAST_ENV"
    printf 'sistr_env\t%s\n' "$SISTR_ENV"
    printf 'amrfinder_env\t%s\n' "$AMRFINDER_ENV"
} > "$PROV"

# ----------------------------------------------------------------------
# Build first-pass rescue candidate set.
#
# Include:
#   - anything that does NOT already have a direct Rigby QC pass
#   - Salmonella fraction >=20%
#
# This includes:
#   - >=70% Salmonella libraries that failed assembly QC
#   - substantial sub-70% mixed libraries
#
# It deliberately does NOT yet define <20% as irrecoverable.
# Those can be assessed later by absolute Salmonella depth.
# ----------------------------------------------------------------------

python3 - "$QC_TABLE" "$CANDIDATES" "$RESCUE_MIN_SALMONELLA_PCT" <<'PY'
import csv
import sys

src, out, minimum = sys.argv[1], sys.argv[2], float(sys.argv[3])

with open(src, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

selected = []

for r in rows:
    pct = float(r["salmonella_read_pct"])

    # Already have a direct, production-quality genome.
    if r["rigby_qc_pass"] == "yes":
        continue

    if pct < minimum:
        continue

    reason = (
        "GATE_PASS_ASSEMBLY_QC_FAIL"
        if r["salmonella_gate_pass"] == "yes"
        else "SUB70_SALMONELLA_RESCUE"
    )

    selected.append({
        "library_id": r["library_id"],
        "biological_isolate_id": r["biological_isolate_id"],
        "sequencing_stream": r["sequencing_stream"],
        "source_type": r["source_type"],
        "specimen_type": r["specimen_type"],
        "dominant_taxon": r["dominant_taxon"],
        "salmonella_read_pct": r["salmonella_read_pct"],
        "original_taxonomy_class": r["taxonomy_class"],
        "original_assembly_length_bp": r["assembly_length_bp"],
        "original_checkm_contamination": r["checkm_contamination"],
        "original_rigby_qc_pass": r["rigby_qc_pass"],
        "rescue_reason": reason,
    })

fields = [
    "library_id",
    "biological_isolate_id",
    "sequencing_stream",
    "source_type",
    "specimen_type",
    "dominant_taxon",
    "salmonella_read_pct",
    "original_taxonomy_class",
    "original_assembly_length_bp",
    "original_checkm_contamination",
    "original_rigby_qc_pass",
    "rescue_reason",
]

with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(selected)

print(f"Rescue candidates: {len(selected)}")

for r in selected:
    print(
        f"{r['library_id']}\t"
        f"{r['salmonella_read_pct']}%\t"
        f"{r['rescue_reason']}"
    )
PY

if [[ "$MODE" == "--plan" ]]; then
    echo
    echo "PLAN ONLY: no rescue analyses launched."
    echo "Candidate table: $CANDIDATES"
    exit 0
fi

process_one() {

    local lib="$1"

    local d="${ROOT}/${lib}"

    local in1="${ANA}/02_host_depleted/${lib}/${lib}_R1.fastq.gz"
    local in2="${ANA}/02_host_depleted/${lib}/${lib}_R2.fastq.gz"

    local report="${d}/${lib}.rescue.report.tsv"

    local kout="${d}/${lib}.kraken.tsv"
    local koutgz="${kout}.gz"

    local all1="${d}/${lib}.salmonella_R1.fastq"
    local all2="${d}/${lib}.salmonella_R2.fastq"

    local norm1="${d}/${lib}.salmonella_120x_R1.fastq.gz"
    local norm2="${d}/${lib}.salmonella_120x_R2.fastq.gz"

    local asmroot="${d}/rescue_assembly"
    local asm="${asmroot}/contigs.fasta"

    local qc="${d}/rescue_qc"
    local td="${d}/rescue_typing"

    local metrics="${d}/rescue_metrics.tsv"

    local kraken_done="${d}/.kraken_complete"
    local extract_done="${d}/.extract_complete"
    local norm_done="${d}/.normalise_complete"
    local lowdepth_flag="${d}/INSUFFICIENT_DEPTH.flag"

    echo
    echo "=============================================================="
    echo "RESCUE: $lib"
    echo "=============================================================="

    [[ -s "$in1" && -s "$in2" ]] || {
        echo "ERROR: Stage-01 reads missing for $lib" >&2
        return 1
    }

    mkdir -p "$d"

    # Existing fully completed rescue — e.g. DEQ172 pilot.
    if [[ \
        -s "$asm" \
        && -s "$qc/quast/report.tsv" \
        && -s "$qc/checkm/checkm_qa.tsv" \
        && -s "$td/sistr.tsv" \
        && -s "$td/mlst.tsv" \
        && -s "$td/amrfinder.tsv" \
    ]]; then
        echo "PASS: existing complete rescued genome + QC + typing found"
        echo "Leaving it untouched."
        return 0
    fi

    if [[ -e "$lowdepth_flag" ]]; then
        echo "PASS: previously classified as insufficient Salmonella depth"
        return 0
    fi

    # --------------------------------------------------------------
    # Kraken2 with per-read classifications.
    # Sentinel prevents an interrupted Kraken run being mistaken
    # for a completed result.
    # --------------------------------------------------------------

    if [[ ! -e "$kraken_done" ]]; then

        rm -f "$report" "$kout" "$koutgz"

        "$KRAKEN2_BIN" \
            --db "$KRAKEN_DB" \
            --threads "$KRAKEN_THREADS" \
            --paired \
            --gzip-compressed \
            --report "$report" \
            --output "$kout" \
            "$in1" "$in2"

        [[ -s "$report" && -s "$kout" ]] || {
            echo "ERROR: incomplete Kraken output for $lib" >&2
            return 1
        }

        touch "$kraken_done"
    fi

    # --------------------------------------------------------------
    # Extract Salmonella taxid 590 + all descendants.
    # --------------------------------------------------------------

    if [[ ! -e "$norm_done" || ! -s "$norm1" || ! -s "$norm2" ]]; then

        if [[ \
            ! -e "$extract_done" \
            || ! -s "$all1" \
            || ! -s "$all2" \
        ]]; then

            rm -f "$all1" "$all2" "$extract_done"

            if [[ ! -s "$kout" && -s "$koutgz" ]]; then
                gzip -dc "$koutgz" > "$kout"
            fi

            [[ -s "$kout" ]] || {
                echo "ERROR: Kraken per-read output missing for $lib" >&2
                return 1
            }

            "$RESCUE_PY" "$KRAKENTOOLS" \
                -k "$kout" \
                -s "$in1" \
                -s2 "$in2" \
                -t 590 \
                -o "$all1" \
                -o2 "$all2" \
                -r "$report" \
                --include-children \
                --fastq-output \
                --noappend

            [[ -s "$all1" && -s "$all2" ]] || {
                echo "ERROR: Salmonella extraction failed for $lib" >&2
                return 1
            }

            touch "$extract_done"
        fi

        # ----------------------------------------------------------
        # Count exact extracted read pairs / bases / target depth.
        # ----------------------------------------------------------

        read -r pairs bases depth < <(
            "$RESCUE_PY" - \
                "$all1" \
                "$all2" \
                "$RESCUE_GENOME_SIZE_BP" <<'PY'
import sys

p1, p2, genome = sys.argv[1], sys.argv[2], float(sys.argv[3])

def count_fastq(path):
    reads = 0
    bases = 0

    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break

            seq = f.readline().rstrip()
            f.readline()
            f.readline()

            reads += 1
            bases += len(seq)

    return reads, bases

n1, b1 = count_fastq(p1)
n2, b2 = count_fastq(p2)

if n1 != n2:
    raise SystemExit(
        f"paired FASTQ count mismatch: R1={n1}, R2={n2}"
    )

total_bases = b1 + b2
depth = total_bases / genome

print(n1, total_bases, depth)
PY
        )

        printf \
            'library_id\textracted_pairs\textracted_bases\traw_salmonella_depth_x\n%s\t%s\t%s\t%.3f\n' \
            "$lib" \
            "$pairs" \
            "$bases" \
            "$depth" \
            > "$metrics"

        # ----------------------------------------------------------
        # Do not assemble a target that has inadequate target depth.
        # This is a rescue-triage rule, not a final genome-QC rule.
        # ----------------------------------------------------------

        if awk \
            -v d="$depth" \
            -v m="$RESCUE_MIN_DEPTH_X" \
            'BEGIN{exit !(d<m)}'
        then
            echo \
                "REVIEW: $lib has only ${depth}x extracted Salmonella depth " \
                "(<${RESCUE_MIN_DEPTH_X}x); no assembly attempted"

            gzip -f "$all1" "$all2"

            if [[ -s "$kout" ]]; then
                gzip -f "$kout"
            fi

            touch "$lowdepth_flag"
            return 0
        fi

        # ----------------------------------------------------------
        # Standardise rescued Salmonella data to <=120x.
        # Same Rasusa pattern used by production.
        # ----------------------------------------------------------

        if awk \
            -v d="$depth" \
            -v c="$RESCUE_CAP_X" \
            'BEGIN{exit !(d>c)}'
        then

            "$RASUSA_BIN" reads \
                --coverage "$RESCUE_CAP_X" \
                --genome-size "$RESCUE_GENOME_SIZE_BP" \
                --seed "$DOWNSAMPLE_SEED" \
                -o "$norm1" \
                -o "$norm2" \
                "$all1" \
                "$all2"

        else

            gzip -c "$all1" > "$norm1"
            gzip -c "$all2" > "$norm2"

        fi

        gzip -t "$norm1" "$norm2"

        touch "$norm_done"

        # Large intermediates are reproducible from Stage-01 reads,
        # Kraken database + report, so do not retain uncompressed FASTQs.
        rm -f "$all1" "$all2" "$extract_done"

        # Retain per-read Kraken classifications compressed for audit.
        if [[ -s "$kout" ]]; then
            gzip -f "$kout"
        fi
    fi

    # --------------------------------------------------------------
    # SPAdes — identical production settings.
    # --------------------------------------------------------------

    if [[ ! -s "$asm" ]]; then

        if [[ -e "$asmroot" ]]; then
            rm -rf "$asmroot"
        fi

        mkdir -p "$asmroot"

        "$SPADES_BIN" \
            --isolate \
            -1 "$norm1" \
            -2 "$norm2" \
            -o "$asmroot/spades" \
            -t "$SPADES_THREADS" \
            -m "$SPADES_MEMORY_GB"

        [[ -s "$asmroot/spades/contigs.fasta" ]] || {
            echo "ERROR: SPAdes produced no contigs for $lib" >&2
            return 1
        }

        cp "$asmroot/spades/contigs.fasta" "$asm"
    fi

    # --------------------------------------------------------------
    # QUAST — identical production settings.
    # --------------------------------------------------------------

    if [[ ! -s "$qc/quast/report.tsv" ]]; then

        rm -rf "$qc/quast"
        mkdir -p "$qc/quast"

        "$CONDA" run --no-capture-output -n "$QUAST_ENV" \
            quast.py \
                --min-contig 500 \
                --threads 2 \
                -o "$qc/quast" \
                "$asm"
    fi

    # --------------------------------------------------------------
    # CheckM — identical production settings.
    # --------------------------------------------------------------

    if [[ ! -s "$qc/checkm/checkm_qa.tsv" ]]; then

        rm -rf "$qc/checkm"
        mkdir -p "$qc/checkm/bins"

        ln -s \
            "$(cd "$(dirname "$asm")" && pwd -P)/$(basename "$asm")" \
            "$qc/checkm/bins/${lib}_RESCUED.fna"

        "$CONDA" run --no-capture-output -n "$CHECKM_ENV" \
            checkm lineage_wf \
                -x fna \
                -t "$CHECKM_THREADS" \
                --pplacer_threads 1 \
                --tab_table \
                -f "$qc/checkm/checkm_qa.tsv" \
                "$qc/checkm/bins" \
                "$qc/checkm/work"
    fi

    # --------------------------------------------------------------
    # Apply CURRENT production final assembly thresholds.
    #
    # No relaxed final QC here.
    # --------------------------------------------------------------

    local qcpass

    qcpass="$(
        python3 - \
            "$qc/quast/report.tsv" \
            "$qc/checkm/checkm_qa.tsv" \
            "$MIN_GENOME_LENGTH" \
            "$MAX_GENOME_LENGTH" \
            "$MAX_CONTIGS" \
            "$MIN_N50" \
            "$MIN_COMPLETENESS" \
            "$MAX_CONTAMINATION" <<'PY'
import csv
import sys

(
    qpath,
    cpath,
    minlen,
    maxlen,
    maxcontigs,
    minn50,
    mincomp,
    maxcontam,
) = sys.argv[1:]

q = {}

with open(qpath) as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) >= 2:
            q[row[0]] = row[1]

with open(cpath, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

if len(rows) != 1:
    raise SystemExit(
        f"expected one CheckM row; found {len(rows)}"
    )

r = rows[0]

def field(name):
    target = name.lower().replace(" ", "")

    for k, v in r.items():
        if k.lower().replace(" ", "") == target:
            return v

    raise KeyError(name)

length = int(q["Total length"])
contigs = int(q["# contigs"])
n50 = int(q["N50"])

completeness = float(field("Completeness"))
contamination = float(field("Contamination"))

passed = (
    length >= int(minlen)
    and length <= int(maxlen)
    and contigs <= int(maxcontigs)
    and n50 >= int(minn50)
    and completeness >= float(mincomp)
    and contamination <= float(maxcontam)
)

print("yes" if passed else "no")
PY
    )"

    if [[ "$qcpass" != "yes" ]]; then
        echo \
            "REVIEW: rescued assembly fails CURRENT production genome QC; " \
            "typing not run"
        return 0
    fi

    # --------------------------------------------------------------
    # SISTR + MLST + AMRFinderPlus
    # identical production methods.
    # --------------------------------------------------------------

    if [[ \
        ! -s "$td/sistr.tsv" \
        || ! -s "$td/mlst.tsv" \
        || ! -s "$td/amrfinder.tsv" \
    ]]; then

        rm -rf "$td"
        mkdir -p "$td"

        "$CONDA" run --no-capture-output -n "$SISTR_ENV" \
            sistr \
                --qc \
                --alleles-output "$td/sistr_alleles.json" \
                --cgmlst-profiles "$td/sistr_cgmlst.csv" \
                -f tab \
                -o "$td/sistr" \
                "$asm"

        [[ -s "$td/sistr.tab" ]] || {
            echo "ERROR: SISTR output missing for $lib" >&2
            return 1
        }

        mv "$td/sistr.tab" "$td/sistr.tsv"

        PATH="$NCBI_BLAST_BIN_DIR:$(dirname "$ANY2FASTA_BIN"):$PATH" \
            "$MLST_BIN" \
            "$asm" \
            > "$td/mlst.tsv"

        [[ -s "$td/mlst.tsv" ]] || {
            echo "ERROR: MLST output missing for $lib" >&2
            return 1
        }

        "$CONDA" run --no-capture-output -n "$AMRFINDER_ENV" \
            amrfinder \
                -n "$asm" \
                -d "$AMRFINDER_DB" \
                --threads 2 \
                --plus \
                -o "$td/amrfinder.tsv"

        [[ -s "$td/amrfinder.tsv" ]] || {
            echo "ERROR: AMRFinderPlus output missing for $lib" >&2
            return 1
        }
    fi

    echo "PASS: rescue processing complete for $lib"
}

while IFS=$'\t' read -r lib _rest
do
    [[ "$lib" == "library_id" ]] && continue

    process_one "$lib"

done < "$CANDIDATES"

# ----------------------------------------------------------------------
# Final consolidated rescue summary.
# ----------------------------------------------------------------------

python3 - \
    "$CANDIDATES" \
    "$ROOT" \
    "$SUMMARY" \
    "$MIN_GENOME_LENGTH" \
    "$MAX_GENOME_LENGTH" \
    "$MAX_CONTIGS" \
    "$MIN_N50" \
    "$MIN_COMPLETENESS" \
    "$MAX_CONTAMINATION" \
    "$RESCUE_MIN_DEPTH_X" <<'PY'
import csv
import sys
from pathlib import Path
from collections import Counter

cand_path, root, out = sys.argv[1:4]

(
    minlen,
    maxlen,
    maxcontigs,
    minn50,
    mincomp,
    maxcontam,
    mindepth,
) = map(float, sys.argv[4:])

root = Path(root)

with open(cand_path, newline="") as f:
    candidates = list(csv.DictReader(f, delimiter="\t"))

def checkm_value(row, name):
    target = name.lower().replace(" ", "")

    for k, v in row.items():
        if k.lower().replace(" ", "") == target:
            return v

    return ""

results = []

for source in candidates:

    lib = source["library_id"]
    d = root / lib

    rec = dict(source)

    rec.update({
        "raw_salmonella_depth_x": "",
        "assembly_length_bp": "",
        "contig_count": "",
        "n50_bp": "",
        "checkm_completeness": "",
        "checkm_contamination": "",
        "assembly_qc_pass": "",
        "sistr_qc_status": "",
        "sistr_serovar": "",
        "sistr_subspecies": "",
        "sistr_cgmlst_st": "",
        "mlst_scheme": "",
        "mlst_st": "",
        "rescue_status": "",
    })

    metrics = d / "rescue_metrics.tsv"

    if metrics.exists():
        with metrics.open(newline="") as f:
            m = list(csv.DictReader(f, delimiter="\t"))

        if m:
            rec["raw_salmonella_depth_x"] = \
                m[0]["raw_salmonella_depth_x"]

    depth = float(rec["raw_salmonella_depth_x"] or 0)

    if depth and depth < mindepth:
        rec["rescue_status"] = "INSUFFICIENT_DEPTH"
        results.append(rec)
        continue

    qpath = d / "rescue_qc/quast/report.tsv"
    cpath = d / "rescue_qc/checkm/checkm_qa.tsv"

    if not (qpath.exists() and cpath.exists()):
        rec["rescue_status"] = "INCOMPLETE"
        results.append(rec)
        continue

    q = {}

    with qpath.open() as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) >= 2:
                q[row[0]] = row[1]

    with cpath.open(newline="") as f:
        checkm_rows = list(csv.DictReader(f, delimiter="\t"))

    if len(checkm_rows) != 1:
        rec["rescue_status"] = "CHECKM_REVIEW"
        results.append(rec)
        continue

    rec["assembly_length_bp"] = q.get("Total length", "")
    rec["contig_count"] = q.get("# contigs", "")
    rec["n50_bp"] = q.get("N50", "")

    rec["checkm_completeness"] = \
        checkm_value(checkm_rows[0], "Completeness")

    rec["checkm_contamination"] = \
        checkm_value(checkm_rows[0], "Contamination")

    assembly_ok = (
        float(rec["assembly_length_bp"]) >= minlen
        and float(rec["assembly_length_bp"]) <= maxlen
        and float(rec["contig_count"]) <= maxcontigs
        and float(rec["n50_bp"]) >= minn50
        and float(rec["checkm_completeness"]) >= mincomp
        and float(rec["checkm_contamination"]) <= maxcontam
    )

    rec["assembly_qc_pass"] = \
        "yes" if assembly_ok else "no"

    if not assembly_ok:
        rec["rescue_status"] = "ASSEMBLY_QC_FAIL"
        results.append(rec)
        continue

    sistr = d / "rescue_typing/sistr.tsv"
    mlst = d / "rescue_typing/mlst.tsv"

    if not (sistr.exists() and mlst.exists()):
        rec["rescue_status"] = "TYPING_INCOMPLETE"
        results.append(rec)
        continue

    with sistr.open(newline="") as f:
        sr = list(csv.DictReader(f, delimiter="\t"))

    if len(sr) == 1:
        rec["sistr_qc_status"] = sr[0].get("qc_status", "")
        rec["sistr_serovar"] = sr[0].get("serovar", "")
        rec["sistr_subspecies"] = \
            sr[0].get("cgmlst_subspecies", "")
        rec["sistr_cgmlst_st"] = \
            sr[0].get("cgmlst_ST", "")

    parts = mlst.read_text().strip().split("\t")

    if len(parts) >= 3:
        rec["mlst_scheme"] = parts[1]
        rec["mlst_st"] = parts[2]

    # A library-level RESCUED_PASS means:
    # - rescued assembly passes unchanged production genome QC
    # - SISTR QC is PASS
    # - MLST identifies Salmonella with an assigned ST
    #
    # Biological-isolate reconciliation remains a later step.
    if (
        rec["sistr_qc_status"] == "PASS"
        and rec["mlst_scheme"] == "salmonella"
        and rec["mlst_st"] not in ("", "-")
    ):
        rec["rescue_status"] = "RESCUED_PASS"
    else:
        rec["rescue_status"] = "TYPING_REVIEW"

    results.append(rec)

if not results:
    raise SystemExit("No rescue candidates found")

fields = list(results[0])

with open(out, "w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
    )

    w.writeheader()
    w.writerows(results)

counts = Counter(
    r["rescue_status"]
    for r in results
)

print()
print("===== RESCUE SUMMARY =====")

for status, n in sorted(counts.items()):
    print(f"{status}: {n}")

print(f"Summary: {out}")
PY
