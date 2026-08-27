#!/opt/homebrew/bin/bash
set -euo pipefail
IFS=$'\n\t'
umask 077

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/tints_config.env"

HELPER="${SCRIPT_DIR}/tints_helpers.py"
MANIFEST="${PROJECT}/00_metadata/TiNTS_library_manifest.tsv"
SOFTWARE_VERSIONS="${PROJECT}/00_metadata/software_versions.tsv"
DATABASE_VERSIONS="${PROJECT}/00_metadata/database_versions.tsv"

ANALYSIS_DIR="${PROJECT}/06_genomics"
RESULTS_DIR="${PROJECT}/RESULTS"
LOG_DIR="${PROJECT}/logs"

MIN_FREE_DISK_GB=200

stamp() {
    printf '\n[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

pass() {
    printf 'PASS: %s\n' "$*"
}

check_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

check_dir() {
    [[ -d "$1" ]] || die "Missing directory: $1"
}

check_executable() {
    [[ -x "$1" ]] || die "Missing/non-executable tool: $1"
}

sha256_of() {
    shasum -a 256 "$1" | awk '{print $1}'
}

require_sha256() {
    local file="$1"
    local expected="$2"
    local observed

    check_file "$file"
    observed="$(sha256_of "$file")"

    [[ "$observed" == "$expected" ]] || \
        die "SHA256 mismatch for $file: expected=$expected observed=$observed"

    pass "SHA256 verified: $file"
}

conda_has_command() {
    local env="$1"
    local command_name="$2"

    "$CONDA" run --no-capture-output -n "$env" \
        bash -lc "command -v '$command_name' >/dev/null" \
        >/dev/null 2>&1 \
        || die "Command '$command_name' unavailable in Conda environment '$env'"

    pass "$env: $command_name"
}

check_database_record() {
    local database="$1"
    local expected_version="$2"
    local expected_hash="$3"

    awk -F'\t' \
        -v db="$database" \
        -v ver="$expected_version" \
        -v hash="$expected_hash" '
        NR > 1 &&
        $1 == db &&
        $2 == ver &&
        $4 == hash {
            found=1
        }
        END {
            exit(found ? 0 : 1)
        }
    ' "$DATABASE_VERSIONS" \
        || die "database_versions.tsv does not contain expected record for $database"

    pass "database_versions.tsv: $database $expected_version"
}

output_has_content() {
    local d="$1"

    [[ -d "$d" ]] || return 1

    find "$d" -mindepth 1 \
        \( -type f -o -type l \) \
        -print -quit 2>/dev/null | grep -q .
}

stamp "TiNTS STAGE 00 PREFLIGHT START"

# ----------------------------------------------------------------------
# 00A. Project/configuration
# ----------------------------------------------------------------------

stamp "00A: project and configuration"

[[ "$PROJECT" == "${HOME}/TiNTS_genomics" ]] \
    || die "Unexpected PROJECT root: $PROJECT"

[[ "$(pwd -P)" == "$(cd "$PROJECT" && pwd -P)" ]] \
    || die "Run this workflow from the TiNTS project root: $PROJECT"

check_file "${SCRIPT_DIR}/tints_config.env"
check_file "$HELPER"

if grep -qE 'LEGACY_ENV|tints-legacy|standard16_20260626|amrfinder/latest' \
    "${SCRIPT_DIR}/tints_config.env"; then
    die "Obsolete monolithic-environment/database reference remains in config"
fi

pass "project root and multi-environment configuration"

# ----------------------------------------------------------------------
# 00B. Frozen metadata
# ----------------------------------------------------------------------

stamp "00B: frozen metadata"

require_sha256 "$LIBRARY_METADATA" "$LIBRARY_METADATA_SHA256"
require_sha256 "$ISOLATE_METADATA" "$ISOLATE_METADATA_SHA256"
require_sha256 "$FULL_METADATA" "$FULL_METADATA_SHA256"

python3 - "$LIBRARY_METADATA" "$ISOLATE_METADATA" "$FULL_METADATA" \
    "$EXPECTED_LIBRARIES" "$EXPECTED_BIOLOGICAL_ISOLATES" \
    "$EXPECTED_MLW_LIBRARIES" "$EXPECTED_CGR_LIBRARIES" \
    "$EXPECTED_HUMAN_ISOLATES" "$EXPECTED_ENVIRONMENTAL_ISOLATES" <<'PY'
import csv
import sys
from collections import Counter

(
    library_path,
    isolate_path,
    full_path,
    expected_libraries,
    expected_isolates,
    expected_mlw,
    expected_cgr,
    expected_human,
    expected_environment,
) = sys.argv[1:]

expected_libraries = int(expected_libraries)
expected_isolates = int(expected_isolates)
expected_mlw = int(expected_mlw)
expected_cgr = int(expected_cgr)
expected_human = int(expected_human)
expected_environment = int(expected_environment)

def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

libraries = read_tsv(library_path)
isolates = read_tsv(isolate_path)
full = read_tsv(full_path)

required_library_columns = {
    "library_id",
    "biological_isolate_id",
    "sequencing_stream",
    "household_id",
    "participant_id",
    "collection_date",
    "source_type",
    "specimen_type",
}

required_isolate_columns = {
    "biological_isolate_id",
    "household_id",
    "participant_id",
    "collection_date",
    "source_type",
    "specimen_type",
}

if not libraries:
    raise SystemExit("library_metadata.tsv contains no data rows")

if not isolates:
    raise SystemExit("isolate_metadata.tsv contains no data rows")

if not required_library_columns.issubset(libraries[0]):
    raise SystemExit("library_metadata.tsv has unexpected/missing columns")

if not required_isolate_columns.issubset(isolates[0]):
    raise SystemExit("isolate_metadata.tsv has unexpected/missing columns")

library_ids = [x["library_id"] for x in libraries]
bio_ids_library = [x["biological_isolate_id"] for x in libraries]
bio_ids_isolate = [x["biological_isolate_id"] for x in isolates]

if len(libraries) != expected_libraries:
    raise SystemExit(
        f"Expected {expected_libraries} library rows; found {len(libraries)}"
    )

if len(set(library_ids)) != expected_libraries:
    raise SystemExit("library_id values are not unique")

if len(isolates) != expected_isolates:
    raise SystemExit(
        f"Expected {expected_isolates} isolate rows; found {len(isolates)}"
    )

if len(set(bio_ids_isolate)) != expected_isolates:
    raise SystemExit("biological_isolate_id values are not unique in isolate metadata")

if len(set(bio_ids_library)) != expected_isolates:
    raise SystemExit(
        "library_metadata.tsv does not resolve to the expected number "
        "of biological isolates"
    )

if set(bio_ids_library) != set(bio_ids_isolate):
    missing = sorted(set(bio_ids_isolate) - set(bio_ids_library))
    extra = sorted(set(bio_ids_library) - set(bio_ids_isolate))
    raise SystemExit(
        f"Biological-isolate mismatch between metadata files. "
        f"Missing={missing}; extra={extra}"
    )

streams = Counter(x["sequencing_stream"] for x in libraries)

if streams != Counter({"MLW": expected_mlw, "CGRdeep": expected_cgr}):
    raise SystemExit(f"Unexpected sequencing-stream counts: {dict(streams)}")

sources = Counter(x["source_type"] for x in isolates)

if sources != Counter(
    {"human": expected_human, "environment": expected_environment}
):
    raise SystemExit(f"Unexpected biological source counts: {dict(sources)}")

if len(full) != expected_libraries:
    raise SystemExit(
        f"Expected {expected_libraries} full_genomics_metadata rows; found {len(full)}"
    )

if "library_id" not in full[0]:
    raise SystemExit("full_genomics_metadata.tsv lacks library_id")

full_ids = [x["library_id"] for x in full]

if len(set(full_ids)) != expected_libraries:
    raise SystemExit("full_genomics_metadata.tsv library_id values are not unique")

if set(full_ids) != set(library_ids):
    missing = sorted(set(library_ids) - set(full_ids))
    extra = sorted(set(full_ids) - set(library_ids))
    raise SystemExit(
        f"Library mismatch between frozen metadata files. "
        f"Missing={missing}; extra={extra}"
    )

print(
    f"PASS: metadata counts: "
    f"libraries={len(libraries)} "
    f"isolates={len(isolates)} "
    f"MLW={streams['MLW']} "
    f"CGRdeep={streams['CGRdeep']} "
    f"human={sources['human']} "
    f"environment={sources['environment']}"
)
PY

# ----------------------------------------------------------------------
# 00C. FASTQ identity and pairing
# ----------------------------------------------------------------------

stamp "00C: metadata-to-FASTQ identity and pairing"

check_dir "$MLW_RAW_DIR"
check_dir "$CGR_ANALYSIS_FASTQ_DIR"

python3 - "$LIBRARY_METADATA" \
    "$MLW_RAW_DIR" \
    "$CGR_ANALYSIS_FASTQ_DIR" \
    "$MANIFEST" <<'PY00C'
import csv
import os
import re
import sys
from pathlib import Path

metadata_path, mlw_raw_dir, cgr_dir, output_path = sys.argv[1:]

with open(metadata_path, newline="") as handle:
    metadata = list(csv.DictReader(handle, delimiter="\t"))

expected = {row["library_id"]: row for row in metadata}

if len(expected) != len(metadata):
    raise SystemExit("Duplicate library_id in library metadata")

observed = {}

# ------------------------------------------------------------
# MLW: derive biological/library ID directly from authoritative
# raw filename. Do not require or create analysis symlinks.
#
# Example:
# 139-DFF12V_B1_GTTGGACGGT-AGGACACTGT_L008_R1.fastq.gz
# -> DFF12V_B1_MLW
# ------------------------------------------------------------

mlw_dir = Path(mlw_raw_dir)

def parse_mlw_r1(filename):
    stem = re.sub(r"\.fastq\.gz$", "", filename, flags=re.I)
    stem = re.sub(r"_R1$", "", stem, flags=re.I)

    m = re.match(r"^(\d+)-(.+)$", stem)
    if not m:
        raise SystemExit(
            f"MLW filename lacks expected numeric-prefix format: {filename}"
        )

    sample_number, remainder = m.groups()

    remainder = re.sub(
        r"_[ACGTN]+(?:-[ACGTN]+)?_L\d{3}$",
        "",
        remainder,
        flags=re.I,
    )
    remainder = re.sub(r"_L\d{3}$", "", remainder, flags=re.I)

    biological_id = remainder.replace("-", "_")
    library_id = f"{biological_id}_MLW"

    return sample_number, biological_id, library_id

mlw_count = 0

for r1 in sorted(mlw_dir.glob("*_R1.fastq.gz")):
    sample_number, biological_id, library_id = parse_mlw_r1(r1.name)

    r2_name = re.sub(
        r"_R1\.fastq\.gz$",
        "_R2.fastq.gz",
        r1.name,
        flags=re.I,
    )
    r2 = mlw_dir / r2_name

    if not r2.is_file():
        raise SystemExit(f"Missing MLW R2 mate for {r1.name}: {r2}")

    if library_id in observed:
        raise SystemExit(
            f"Duplicate parsed MLW library ID {library_id}; "
            f"filename parsing is not one-to-one"
        )

    observed[library_id] = {
        "sequencing_stream": "MLW",
        "biological_isolate_id": biological_id,
        "r1": str(r1.resolve()),
        "r2": str(r2.resolve()),
        "original_sample_number": sample_number,
    }
    mlw_count += 1

# Catch orphan MLW R2 files independently.
mlw_r1_names = {
    re.sub(r"_R1\.fastq\.gz$", "", x.name, flags=re.I)
    for x in mlw_dir.glob("*_R1.fastq.gz")
}
mlw_r2_names = {
    re.sub(r"_R2\.fastq\.gz$", "", x.name, flags=re.I)
    for x in mlw_dir.glob("*_R2.fastq.gz")
}

if mlw_r1_names != mlw_r2_names:
    missing_r2 = sorted(mlw_r1_names - mlw_r2_names)
    missing_r1 = sorted(mlw_r2_names - mlw_r1_names)
    raise SystemExit(
        f"MLW pairing mismatch. Missing_R2={missing_r2}; "
        f"Missing_R1={missing_r1}"
    )

# ------------------------------------------------------------
# CGRdeep: converted analysis FASTQs already use canonical
# library IDs, e.g. <biological_id>_CGRdeep_R1.fastq.gz.
# ------------------------------------------------------------

cgr_path = Path(cgr_dir)
cgr_count = 0

for r1 in sorted(cgr_path.glob("*_R1.fastq.gz")):
    library_id = r1.name[:-len("_R1.fastq.gz")]
    r2 = cgr_path / f"{library_id}_R2.fastq.gz"

    if not r2.is_file():
        raise SystemExit(f"Missing CGRdeep R2 mate for {library_id}: {r2}")

    if library_id in observed:
        raise SystemExit(f"Duplicate observed library ID: {library_id}")

    observed[library_id] = {
        "sequencing_stream": "CGRdeep",
        "biological_isolate_id": (
            library_id[:-len("_CGRdeep")]
            if library_id.endswith("_CGRdeep")
            else library_id
        ),
        "r1": str(r1.resolve()),
        "r2": str(r2.resolve()),
        "original_sample_number": "NA",
    }
    cgr_count += 1

cgr_r1_names = {
    x.name[:-len("_R1.fastq.gz")]
    for x in cgr_path.glob("*_R1.fastq.gz")
}
cgr_r2_names = {
    x.name[:-len("_R2.fastq.gz")]
    for x in cgr_path.glob("*_R2.fastq.gz")
}

if cgr_r1_names != cgr_r2_names:
    missing_r2 = sorted(cgr_r1_names - cgr_r2_names)
    missing_r1 = sorted(cgr_r2_names - cgr_r1_names)
    raise SystemExit(
        f"CGRdeep pairing mismatch. Missing_R2={missing_r2}; "
        f"Missing_R1={missing_r1}"
    )

# ------------------------------------------------------------
# Frozen metadata is authoritative.
# ------------------------------------------------------------

expected_ids = set(expected)
observed_ids = set(observed)

missing = sorted(expected_ids - observed_ids)
unexpected = sorted(observed_ids - expected_ids)

if missing or unexpected:
    raise SystemExit(
        "Metadata/FASTQ library identity mismatch. "
        f"Missing={missing}; unexpected={unexpected}"
    )

rows = []

for library_id in sorted(expected):
    meta = expected[library_id]
    obs = observed[library_id]

    if meta["sequencing_stream"] != obs["sequencing_stream"]:
        raise SystemExit(
            f"Sequencing-stream mismatch for {library_id}: "
            f"metadata={meta['sequencing_stream']} "
            f"FASTQ={obs['sequencing_stream']}"
        )

    if meta["biological_isolate_id"] != obs["biological_isolate_id"]:
        raise SystemExit(
            f"Biological-isolate mismatch for {library_id}: "
            f"metadata={meta['biological_isolate_id']} "
            f"FASTQ-derived={obs['biological_isolate_id']}"
        )

    if os.path.getsize(obs["r1"]) == 0 or os.path.getsize(obs["r2"]) == 0:
        raise SystemExit(f"Empty FASTQ file for {library_id}")

    rows.append({
        **meta,
        "r1": obs["r1"],
        "r2": obs["r2"],
    })

fields = list(metadata[0].keys()) + ["r1", "r2"]

Path(output_path).parent.mkdir(parents=True, exist_ok=True)
tmp = output_path + ".tmp"

with open(tmp, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

os.replace(tmp, output_path)

print(
    f"PASS: FASTQ discovery: MLW={mlw_count} CGRdeep={cgr_count} "
    f"total={len(observed)}"
)
print(
    f"PASS: expected metadata libraries={len(expected)}; "
    f"observed paired FASTQ libraries={len(observed)}; "
    f"missing=0; unexpected=0"
)
print(f"PASS: deterministic manifest written: {output_path}")
PY00C

manifest_rows="$(awk 'END{print NR-1}' "$MANIFEST")"
[[ "$manifest_rows" -eq "$EXPECTED_LIBRARIES" ]] \
    || die "Manifest row count is $manifest_rows, expected $EXPECTED_LIBRARIES"

pass "manifest contains $manifest_rows paired libraries"

# ----------------------------------------------------------------------
# 00D. FASTQ gzip integrity
# ----------------------------------------------------------------------

stamp "00D: FASTQ gzip integrity"

python3 - "$MANIFEST" "$EXPECTED_LIBRARIES" <<'PY00D'
import csv
import gzip
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
expected_libraries = int(sys.argv[2])

with manifest_path.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

if len(rows) != expected_libraries:
    raise SystemExit(
        f"Manifest contains {len(rows)} libraries; "
        f"expected {expected_libraries}"
    )

checked_files = 0

for row in rows:
    library_id = row["library_id"]

    for column in ("r1", "r2"):
        value = row.get(column, "")
        if not value:
            raise SystemExit(
                f"Missing {column} path in manifest for {library_id}"
            )

        path = Path(value)

        if not path.is_file():
            raise SystemExit(
                f"FASTQ does not exist for {library_id} {column}: {path}"
            )

        # Read the entire compressed stream so CRC/truncation errors are detected.
        try:
            with gzip.open(path, "rb") as handle:
                while handle.read(8 * 1024 * 1024):
                    pass
        except Exception as exc:
            raise SystemExit(
                f"gzip integrity failure for {library_id} {column}: "
                f"{path}: {exc}"
            )

        checked_files += 1

print(
    f"PASS: gzip integrity verified for "
    f"{checked_files} FASTQ files across {len(rows)} libraries"
)
PY00D

# ----------------------------------------------------------------------
# 00E. Native/source software
# ----------------------------------------------------------------------

stamp "00E: native and source-installed software"

check_executable "$FASTP_BIN"
check_executable "$BOWTIE2_BIN"
check_executable "$KRAKEN2_BIN"
check_executable "$SPADES_BIN"
check_executable "$SEQKIT_BIN"
check_executable "$SAMTOOLS_BIN"
check_executable "$RASUSA_BIN"
check_executable "$MAFFT_BIN"
check_executable "$MLST_BIN"
check_executable "$ANY2FASTA_BIN"
check_executable "$IQTREE3_BIN"

pass "all required native/source executable paths resolve"

require_sha256 "$IQTREE3_BIN" "$IQTREE3_SHA256"

check_file "$MLST_ENV_FILE"
check_file "$BLAST_ENV_FILE"
check_dir "$NCBI_BLAST_BIN_DIR"

for blast_tool in blastn makeblastdb; do
    check_executable "${NCBI_BLAST_BIN_DIR}/${blast_tool}"
done

pass "pinned NCBI BLAST+ toolchain resolves"

# ----------------------------------------------------------------------
# 00F. Dedicated Conda environments
# ----------------------------------------------------------------------

stamp "00F: dedicated Conda environments"

check_executable "$CONDA"

conda_has_command "$CHECKM_ENV" checkm
conda_has_command "$QUAST_ENV" quast.py
conda_has_command "$SISTR_ENV" sistr
conda_has_command "$AMRFINDER_ENV" amrfinder
conda_has_command "$SNP_TOOLS_ENV" snp-sites
conda_has_command "$SNP_TOOLS_ENV" snp-dists
conda_has_command "$GUBBINS_ENV" run_gubbins.py
conda_has_command "$GUBBINS_ENV" iqtree

"$CONDA" run -n "$PLOT_ENV" python -c 'import matplotlib, Bio; from Bio import Phylo'     >/dev/null 2>&1     || die "Plotting environment lacks matplotlib and/or Bio.Phylo"

pass "$PLOT_ENV: matplotlib + Bio.Phylo"

declare -A ENV_HASHES=(
    ["$CHECKM_ENV_EXPORT"]="3ac3bc3f383f236ee0db8763af649392dea111a934015cfc90eeea0ff10e6351"
    ["$QUAST_ENV_EXPORT"]="35696cbe3d9e889ed990c8f470fc7d41eb4acf7c6e6753250e3963a0ec3b11f1"
    ["$SISTR_ENV_EXPORT"]="1fe84976022ddafa5b4607a38dff58a91fabe2eedeaf2efde59c7933af2f8258"
    ["$AMRFINDER_ENV_EXPORT"]="9dc873666d70b5adb64239a99fc5d5e5f12972dafa0922b6c8e04d64f6fae8c7"
    ["$SNP_TOOLS_ENV_EXPORT"]="a9e2d64c8b756063ff605993bbc7f46df49a512fcb0831e5a212c1fe0e90dc24"
    ["$GUBBINS_ENV_EXPORT"]="5181136dd3d6528bbc12a3aa5a2488e852d53b97affe90645b891d7b09cc6642"
    ["$PLOT_ENV_EXPORT"]="$PLOT_ENV_EXPORT_SHA256"
)

for env_file in "${!ENV_HASHES[@]}"; do
    require_sha256 "$env_file" "${ENV_HASHES[$env_file]}"
done

# ----------------------------------------------------------------------
# 00G. Databases and human reference
# ----------------------------------------------------------------------

stamp "00G: databases and human reference"

check_file "$DATABASE_VERSIONS"

check_database_record \
    "CheckM" \
    "$CHECKM_DB_SNAPSHOT" \
    "$CHECKM_DB_MANIFEST_SHA256"

check_database_record \
    "AMRFinderPlus" \
    "$AMRFINDER_DB_VERSION" \
    "$AMRFINDER_DB_MANIFEST_SHA256"

check_database_record \
    "$HUMAN_REFERENCE_NAME" \
    "$HUMAN_REFERENCE_ACCESSION" \
    "$HUMAN_REFERENCE_MANIFEST_SHA256"

check_database_record \
    "Kraken2 Standard-16" \
    "$KRAKEN_DB_SNAPSHOT" \
    "$KRAKEN_DB_MANIFEST_SHA256"

check_dir "$CHECKM_DB"
check_dir "$AMRFINDER_DB"
check_dir "$KRAKEN_DB"
check_file "$HUMAN_REFERENCE_FASTA"

for f in hash.k2d opts.k2d taxo.k2d; do
    check_file "${KRAKEN_DB}/${f}"
done

bt2_count="$(
    find "$(dirname "$HUMAN_BOWTIE2_INDEX")" \
        -maxdepth 1 \
        -type f \
        -name "$(basename "$HUMAN_BOWTIE2_INDEX").*.bt2*" \
        | wc -l | tr -d ' '
)"

[[ "$bt2_count" -eq 6 ]] \
    || die "Expected 6 Bowtie2 index files; found $bt2_count"

pass "six Bowtie2 human-reference index files present"

# ----------------------------------------------------------------------
# 00H. Docker and pinned containers
# ----------------------------------------------------------------------

stamp "00H: Docker and pinned container images"

command -v docker >/dev/null \
    || die "Docker CLI unavailable"

docker info >/dev/null 2>&1 \
    || die "Docker Desktop/server is not operational"

pass "Docker server operational"

docker image inspect "$PROKKA_IMAGE" >/dev/null 2>&1 \
    || die "Pinned Prokka image not present locally: $PROKKA_IMAGE"

docker image inspect "$PANAROO_IMAGE" >/dev/null 2>&1 \
    || die "Pinned Panaroo image not present locally: $PANAROO_IMAGE"

docker image inspect "$SNIPPY_IMAGE" >/dev/null 2>&1 \
    || die "Pinned Snippy image not present locally: $SNIPPY_IMAGE"

prokka_digests="$(
    docker image inspect \
        --format '{{join .RepoDigests "\n"}}' \
        "$PROKKA_IMAGE" 2>/dev/null || true
)"

panaroo_digests="$(
    docker image inspect \
        --format '{{join .RepoDigests "\n"}}' \
        "$PANAROO_IMAGE" 2>/dev/null || true
)"

snippy_digests="$(
    docker image inspect \
        --format '{{join .RepoDigests "\n"}}' \
        "$SNIPPY_IMAGE" 2>/dev/null || true
)"

grep -Fq "$PROKKA_IMAGE_DIGEST" <<< "$prokka_digests" \
    || die "Local Prokka image digest does not match pinned digest"

grep -Fq "$PANAROO_IMAGE_DIGEST" <<< "$panaroo_digests" \
    || die "Local Panaroo image digest does not match pinned digest"

grep -Fq "$SNIPPY_IMAGE_DIGEST" <<< "$snippy_digests" \
    || die "Local Snippy image digest does not match pinned digest"

pass "pinned Prokka image digest verified"
pass "pinned Panaroo image digest verified"
pass "pinned Snippy image digest verified"

# ----------------------------------------------------------------------
# 00I. Output-collision guard
# ----------------------------------------------------------------------

stamp "00I: production-output collision guard"

RUN_MODE="${1:-stage00}"

case "$RUN_MODE" in
    stage00|"")
        if output_has_content "$ANALYSIS_DIR"; then
            die "Existing production content detected under $ANALYSIS_DIR."
        fi
        if output_has_content "$RESULTS_DIR"; then
            die "Existing production content detected under $RESULTS_DIR."
        fi
        pass "no existing production-file collision detected"
        ;;
    --production)
        if output_has_content "$ANALYSIS_DIR"; then
            die "Fresh production requested but existing content detected under $ANALYSIS_DIR. Use --resume only for a deliberately resumed validated run."
        fi
        if output_has_content "$RESULTS_DIR"; then
            die "Fresh production requested but existing content detected under $RESULTS_DIR."
        fi
        pass "fresh-production output guard"
        ;;
    --resume)
        if output_has_content "$RESULTS_DIR"; then
            die "RESULTS already contains files; final-results resume is not yet implemented."
        fi
        pass "resume mode: existing analysis output permitted and will require per-stage validation"
        ;;
    *)
        die "Unknown argument: $RUN_MODE. Allowed: no argument, --production, --resume"
        ;;
esac

# ----------------------------------------------------------------------
# 00J. Available disk space
# ----------------------------------------------------------------------

stamp "00J: disk-space check"

free_gb="$(
    df -Pk "$PROJECT" |
    awk 'NR==2 {printf "%d", $4/1024/1024}'
)"

[[ "$free_gb" =~ ^[0-9]+$ ]] \
    || die "Could not determine available disk space"

printf 'Available project-volume space: %s GiB\n' "$free_gb"
printf 'Operational preflight floor: %s GiB\n' "$MIN_FREE_DISK_GB"

(( free_gb >= MIN_FREE_DISK_GB )) \
    || die "Available disk space is below the ${MIN_FREE_DISK_GB} GiB preflight floor"

pass "disk-space guard"

# ----------------------------------------------------------------------
# 00K. Provenance capture
# ----------------------------------------------------------------------

stamp "00K: pre-analysis provenance capture"

check_file "$SOFTWARE_VERSIONS"
check_file "$DATABASE_VERSIONS"

INPUT_CHECKSUMS="${PROJECT}/00_metadata/input_checksums.sha256"

python3 -     "$PROJECT"     "$MLW_RAW_DIR"     "$CGR_RAW_CRAM_DIR"     "$CGR_ANALYSIS_FASTQ_DIR"     "$INPUT_CHECKSUMS" <<'PY00CHECKSUM'
import hashlib
import os
import sys
from pathlib import Path

project, mlw_dir, cgr_cram_dir, cgr_fastq_dir, output = sys.argv[1:]

project = Path(project).resolve()

files = []
files.extend(Path(mlw_dir).glob("*.fastq.gz"))
files.extend(Path(cgr_cram_dir).glob("*.cram"))
files.extend(Path(cgr_fastq_dir).glob("*.fastq.gz"))

files = sorted({x.resolve() for x in files})

expected = 142 + 120 + 80

if len(files) != expected:
    raise SystemExit(
        f"Expected {expected} provenance input files "
        f"(142 MLW FASTQ + 120 CGR CRAM + 80 CGR FASTQ); found {len(files)}"
    )

tmp = Path(str(output) + ".tmp")

with tmp.open("w") as out:
    for path in files:
        h = hashlib.sha256()

        with path.open("rb") as f:
            for block in iter(lambda: f.read(16 * 1024 * 1024), b""):
                h.update(block)

        try:
            rel = path.relative_to(project)
            display = str(rel)
        except ValueError:
            display = str(path)

        out.write(f"{h.hexdigest()}  {display}\n")

os.replace(tmp, output)

print(
    f"PASS: SHA256 provenance captured for {len(files)} "
    f"source/analysis-input files"
)
PY00CHECKSUM

check_file "$INPUT_CHECKSUMS"

PROVENANCE_DIR="$(
    mktemp -d "${PROJECT}/00_metadata/stage00_provenance_$(date +%Y%m%d_%H%M%S)_XXXXXX"
)"

cp "$SOFTWARE_VERSIONS" "$PROVENANCE_DIR/software_versions.tsv"
cp "$DATABASE_VERSIONS" "$PROVENANCE_DIR/database_versions.tsv"

cp "$CHECKM_ENV_EXPORT" "$PROVENANCE_DIR/"
cp "$QUAST_ENV_EXPORT" "$PROVENANCE_DIR/"
cp "$SISTR_ENV_EXPORT" "$PROVENANCE_DIR/"
cp "$AMRFINDER_ENV_EXPORT" "$PROVENANCE_DIR/"
cp "$SNP_TOOLS_ENV_EXPORT" "$PROVENANCE_DIR/"
cp "$GUBBINS_ENV_EXPORT" "$PROVENANCE_DIR/"
cp "$PLOT_ENV_EXPORT" "$PROVENANCE_DIR/"
cp "$MLST_ENV_FILE" "$PROVENANCE_DIR/"
cp "$BLAST_ENV_FILE" "$PROVENANCE_DIR/"

{
    shasum -a 256 \
        "$LIBRARY_METADATA" \
        "$ISOLATE_METADATA" \
        "$FULL_METADATA" \
        "$MANIFEST" \
        "${SCRIPT_DIR}/tints_config.env" \
        "${BASH_SOURCE[0]}" \
        "$HELPER"
} > "$PROVENANCE_DIR/project_file_checksums.sha256"

{
    echo "generated=$(date -Iseconds)"
    echo "project=$PROJECT"
    echo "hostname=$(hostname)"
    echo "machine=$(uname -m)"
    echo "macos=$(sw_vers -productVersion)"
    echo "free_disk_gb=$free_gb"
    echo "expected_libraries=$EXPECTED_LIBRARIES"
    echo "expected_biological_isolates=$EXPECTED_BIOLOGICAL_ISOLATES"
    echo "expected_mlw_libraries=$EXPECTED_MLW_LIBRARIES"
    echo "expected_cgrdeep_libraries=$EXPECTED_CGR_LIBRARIES"
    echo "expected_human_isolates=$EXPECTED_HUMAN_ISOLATES"
    echo "expected_environmental_isolates=$EXPECTED_ENVIRONMENTAL_ISOLATES"
    echo "salmonella_gate_pct=$SALMONELLA_READ_PCT_THRESHOLD"
    echo "depth_cap_x=$DEPTH_CAP_X"
    echo "planning_genome_size_bp=$PLANNING_GENOME_SIZE_BP"
    echo "downsample_seed=$DOWNSAMPLE_SEED"
    echo "human_reference=$HUMAN_REFERENCE_NAME"
    echo "human_reference_accession=$HUMAN_REFERENCE_ACCESSION"
    echo "host_filter_preset=$HOST_FILTER_PRESET"
    echo "kraken_memory_mapping_required=no"
} > "$PROVENANCE_DIR/stage00_parameters.txt"

docker info \
    --format 'server={{.ServerVersion}} os={{.OperatingSystem}} arch={{.Architecture}}' \
    > "$PROVENANCE_DIR/docker_runtime.txt"

"$CONDA" env list > "$PROVENANCE_DIR/conda_env_list.txt"

pass "provenance captured: $PROVENANCE_DIR"

# ----------------------------------------------------------------------
# STOP BY DESIGN
# ----------------------------------------------------------------------

stamp "STAGE 00 COMPLETE"

echo
echo "=============================================================="
echo "STAGE 00 PASSED."
echo "Stage 00 validation completed; downstream analysis stages were NOT launched in this run."
echo "No fastp/Bowtie2/Kraken2/SPAdes analysis has been launched."
echo "Manifest:   $MANIFEST"
echo "Provenance: $PROVENANCE_DIR"
echo "=============================================================="

# ======================================================================
# PRODUCTION CONTINUATION GUARD
# ======================================================================

if [[ "${RUN_MODE:-stage00}" == "stage00" || -z "${RUN_MODE:-}" ]]; then
    echo
    echo "STOPPED AFTER STAGE 00 BY DEFAULT."
    echo "Use --production for a fresh run or --resume for a validated partial run."
    exit 0
fi

[[ "$RUN_MODE" == "--production" || "$RUN_MODE" == "--resume" ]]     || die "Invalid downstream run mode: $RUN_MODE"

ANA="$ANALYSIS_DIR"
STATE="${ANA}/00_state"

mkdir -p \
    "$STATE" \
    "${ANA}/01_read_qc" \
    "${ANA}/02_host_depleted" \
    "${ANA}/03_taxonomy" \
    "${ANA}/04_normalised"

CURRENT_TMP=""

cleanup_tmp() {
    if [[ -n "${CURRENT_TMP:-}" && -d "${CURRENT_TMP:-}" ]]; then
        rm -rf "$CURRENT_TMP"
    fi
}

trap cleanup_tmp EXIT INT TERM

# ======================================================================
# STAGE 01: FASTP + HOST DEPLETION
# ======================================================================

stamp "STAGE 01: fastp + T2T-CHM13v2.0 host depletion"

while IFS=$'\t' read -r lib r1 r2
do
    stamp "Stage 01: $lib"

    qd="${ANA}/01_read_qc/${lib}"
    hd="${ANA}/02_host_depleted/${lib}"

    mkdir -p "$qd" "$hd"
    chmod 700 "$qd" "$hd" 2>/dev/null || true

    out1="${hd}/${lib}_R1.fastq.gz"
    out2="${hd}/${lib}_R2.fastq.gz"
    metrics="${hd}/host_depletion_metrics.tsv"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        fastp_json="${qd}/${lib}.json"

        if [[ -s "$out1" && -s "$out2" && -s "$metrics" && -s "$fastp_json" ]]; then
            gzip -t "$out1" "$out2" \
                || die "Resume validation: corrupt Stage 01 FASTQ for $lib"

            n1="$("$SEQKIT_BIN" stats -T "$out1" | awk -F'\t' 'NR==2{gsub(/,/,"",$4);print $4}')"
            n2="$("$SEQKIT_BIN" stats -T "$out2" | awk -F'\t' 'NR==2{gsub(/,/,"",$4);print $4}')"

            [[ "$n1" =~ ^[0-9]+$ && "$n2" =~ ^[0-9]+$ && "$n1" -eq "$n2" ]] \
                || die "Resume validation: Stage 01 pairing invalid for $lib"

            python3 - "$metrics" "$lib" <<'PYRESUME01'
import csv
import sys

path, expected_lib = sys.argv[1:]

with open(path, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

if len(rows) != 1:
    raise SystemExit(
        f"Resume validation: expected one host-depletion metrics row for {expected_lib}"
    )

r = rows[0]

if r.get("library_id") != expected_lib:
    raise SystemExit(
        f"Resume validation: metrics library mismatch for {expected_lib}"
    )

if r.get("pair_retention_rule") != "BOTH_MATES_UNMAPPED_REQUIRED":
    raise SystemExit(
        f"Resume validation: wrong pair-retention rule for {expected_lib}"
    )
PYRESUME01

            pass "resume: validated existing Stage 01 output for $lib"
            continue
        fi

        if [[ -e "$out1" || -e "$out2" || -e "$metrics" || -e "$fastp_json" ]]; then
            die "Resume validation: partial/incomplete Stage 01 output for $lib"
        fi

    else
        [[ ! -e "$out1" && ! -e "$out2" && ! -e "$metrics" ]] \
            || die "Unexpected pre-existing Stage 01 output for $lib"
    fi

    CURRENT_TMP="$(mktemp -d "${STATE}/prehost_${lib}.XXXXXX")"

    pre1="${CURRENT_TMP}/${lib}_prehost_R1.fastq.gz"
    pre2="${CURRENT_TMP}/${lib}_prehost_R2.fastq.gz"

    "$FASTP_BIN" \
        -i "$r1" \
        -I "$r2" \
        -o "$pre1" \
        -O "$pre2" \
        --detect_adapter_for_pe \
        --thread "$FASTP_THREADS" \
        --json "${qd}/${lib}.json" \
        --html "${qd}/${lib}.html" \
        > "${qd}/fastp.stdout.log" \
        2> "${qd}/fastp.stderr.log"

    gzip -t "$pre1" "$pre2"

    "$BOWTIE2_BIN" "$HOST_FILTER_PRESET" \
        -x "$HUMAN_BOWTIE2_INDEX" \
        -1 "$pre1" \
        -2 "$pre2" \
        -p "$HOST_FILTER_THREADS" \
        2> "${hd}/bowtie2_human.stderr.log" \
    | "$SAMTOOLS_BIN" view \
        -@ "$SAMTOOLS_THREADS" \
        -u \
        -f 12 \
        -F 2304 \
        - \
    | "$SAMTOOLS_BIN" collate \
        -@ "$SAMTOOLS_THREADS" \
        -u \
        -O \
        - \
    | "$SAMTOOLS_BIN" fastq \
        -@ "$SAMTOOLS_THREADS" \
        -1 "$out1" \
        -2 "$out2" \
        -0 /dev/null \
        -s /dev/null \
        -n \
        - \
        >/dev/null

    gzip -t "$out1" "$out2"

    IFS=$'\t' read -r input_reads input_bases < <(
        python3 - "${qd}/${lib}.json" <<'PYFASTP'
import json
import sys

with open(sys.argv[1]) as f:
    j = json.load(f)

s = j["summary"]["after_filtering"]

print(
    int(s["total_reads"]),
    int(s["total_bases"]),
    sep="\t",
)
PYFASTP
    )

    [[ "$input_reads" =~ ^[0-9]+$ ]] \
        || die "Could not parse fastp read count for $lib"

    [[ "$input_bases" =~ ^[0-9]+$ ]] \
        || die "Could not parse fastp base count for $lib"

    (( input_reads % 2 == 0 )) \
        || die "Odd total read count after fastp for $lib"

    input_pairs=$(( input_reads / 2 ))

    IFS=$'\t' read -r r1_n r1_b < <(
        "$SEQKIT_BIN" stats -T "$out1" |
        awk -F'\t' '
            NR==2 {
                gsub(/,/,"",$4)
                gsub(/,/,"",$5)
                print $4 "\t" $5
            }
        '
    )

    IFS=$'\t' read -r r2_n r2_b < <(
        "$SEQKIT_BIN" stats -T "$out2" |
        awk -F'\t' '
            NR==2 {
                gsub(/,/,"",$4)
                gsub(/,/,"",$5)
                print $4 "\t" $5
            }
        '
    )

    [[ "$r1_n" =~ ^[0-9]+$ && "$r2_n" =~ ^[0-9]+$ ]] \
        || die "Could not count host-depleted reads for $lib"

    [[ "$r1_n" -eq "$r2_n" ]] \
        || die "Host depletion broke pairing for $lib: R1=$r1_n R2=$r2_n"

    retained_pairs="$r1_n"
    retained_bases=$(( r1_b + r2_b ))
    removed_pairs=$(( input_pairs - retained_pairs ))

    (( removed_pairs >= 0 )) \
        || die "Negative removed-pair count for $lib"

    removed_pct="$(
        awk -v x="$removed_pairs" -v n="$input_pairs" \
            'BEGIN{if(n==0)print "0.000000";else printf "%.6f",100*x/n}'
    )"

    retained_pct="$(
        awk -v x="$retained_pairs" -v n="$input_pairs" \
            'BEGIN{if(n==0)print "0.000000";else printf "%.6f",100*x/n}'
    )"

    {
        printf '%s\n' \
'library_id	input_reads_after_fastp	input_pairs_after_fastp	input_bases_after_fastp	retained_nonhuman_pairs	removed_pairs_any_human_alignment	removed_pair_pct	retained_pair_pct	retained_nonhuman_bases	bowtie2_preset	human_reference	human_reference_accession	pair_retention_rule'

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$lib" \
            "$input_reads" \
            "$input_pairs" \
            "$input_bases" \
            "$retained_pairs" \
            "$removed_pairs" \
            "$removed_pct" \
            "$retained_pct" \
            "$retained_bases" \
            "$HOST_FILTER_PRESET" \
            "$HUMAN_REFERENCE_NAME" \
            "$HUMAN_REFERENCE_ACCESSION" \
            'BOTH_MATES_UNMAPPED_REQUIRED'
    } > "$metrics"

    rm -rf "$CURRENT_TMP"
    CURRENT_TMP=""

done < <(
    python3 -c 'import csv,sys; f=open(sys.argv[1],newline=""); [print("\t".join((r["library_id"],r["r1"],r["r2"]))) for r in csv.DictReader(f,delimiter="\t")]' "$MANIFEST"
)

python3 "$HELPER" host-summary \
    --manifest "$MANIFEST" \
    --host "${ANA}/02_host_depleted" \
    --genome-size "$PLANNING_GENOME_SIZE_BP" \
    --cap "$DEPTH_CAP_X" \
    --output "${ANA}/02_host_depleted/host_depletion_summary.tsv"

stamp "STAGE 01 COMPLETE"

# ======================================================================
# STAGE 02: KRAKEN2 ON FULL HOST-DEPLETED READS
# ======================================================================

stamp "STAGE 02: Kraken2 taxonomy before any depth normalisation"

while IFS=$'\t' read -r lib
do
    stamp "Stage 02: $lib"

    d="${ANA}/03_taxonomy/${lib}"
    mkdir -p "$d"

    report="${d}/${lib}.report.tsv"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        if [[ -s "$report" ]]; then
            python3 - "$report" "$lib" <<'PYRESUME02'
import sys

path, lib = sys.argv[1:]
have_root = False
have_unclassified = False

with open(path) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue

        rank = p[3].strip()
        taxid = p[4].strip()

        if rank == "R" and taxid == "1":
            have_root = True

        if rank == "U" and taxid == "0":
            have_unclassified = True

if not have_root or not have_unclassified:
    raise SystemExit(
        f"Resume validation: incomplete Kraken2 report for {lib}"
    )
PYRESUME02

            pass "resume: validated existing Stage 02 report for $lib"
            continue
        fi

        if [[ -e "$report" ]]; then
            die "Resume validation: empty/incomplete Stage 02 report for $lib"
        fi

    else
        [[ ! -e "$report" ]] \
            || die "Unexpected pre-existing Kraken2 report for $lib"
    fi

    in1="${ANA}/02_host_depleted/${lib}/${lib}_R1.fastq.gz"
    in2="${ANA}/02_host_depleted/${lib}/${lib}_R2.fastq.gz"

    [[ -s "$in1" && -s "$in2" ]] \
        || die "Missing Stage 01 reads for $lib"

    "$KRAKEN2_BIN" \
        --db "$KRAKEN_DB" \
        --threads "$KRAKEN_THREADS" \
        --paired \
        --gzip-compressed \
        --report "$report" \
        --output /dev/null \
        "$in1" "$in2"

done < <(
    python3 -c 'import csv,sys; f=open(sys.argv[1],newline=""); [print(r["library_id"]) for r in csv.DictReader(f,delimiter="\t")]' "$MANIFEST"
)

python3 "$HELPER" kraken-summary \
    --manifest "$MANIFEST" \
    --kraken "${ANA}/03_taxonomy" \
    --salmonella-taxid "$SALMONELLA_TAXID" \
    --human-taxid "$HUMAN_TAXID" \
    --threshold "$SALMONELLA_READ_PCT_THRESHOLD" \
    --output "${ANA}/03_taxonomy/taxonomy_summary.tsv"

stamp "STAGE 02 COMPLETE"

# ======================================================================
# STAGE 03: 120x DEPTH CAP
# ======================================================================

stamp "STAGE 03: deterministic post-taxonomy depth normalisation"

PLAN="${ANA}/02_host_depleted/host_depletion_summary.tsv"

while IFS=$'\t' read -r lib action
do
    stamp "Stage 03: $lib [$action]"

    d="${ANA}/04_normalised/${lib}"
    mkdir -p "$d"

    in1="${ANA}/02_host_depleted/${lib}/${lib}_R1.fastq.gz"
    in2="${ANA}/02_host_depleted/${lib}/${lib}_R2.fastq.gz"

    out1="${d}/${lib}_R1.fastq.gz"
    out2="${d}/${lib}_R2.fastq.gz"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        if [[ -e "$out1" && -e "$out2" ]]; then
            gzip -t "$out1" "$out2" \
                || die "Resume validation: corrupt Stage 03 FASTQ for $lib"

            n1="$("$SEQKIT_BIN" stats -T "$out1" | awk -F'\t' 'NR==2{gsub(/,/,"",$4);print $4}')"
            n2="$("$SEQKIT_BIN" stats -T "$out2" | awk -F'\t' 'NR==2{gsub(/,/,"",$4);print $4}')"

            [[ "$n1" =~ ^[0-9]+$ && "$n2" =~ ^[0-9]+$ && "$n1" -eq "$n2" ]] \
                || die "Resume validation: Stage 03 pairing invalid for $lib"

            if [[ "$action" == "unchanged" ]]; then
                [[ -L "$out1" && -L "$out2" ]] \
                    || die "Resume validation: unchanged Stage 03 outputs are not symlinks for $lib"

                resolved_out1="$(cd "$(dirname "$out1")" && cd "$(dirname "$(readlink "$out1")")" 2>/dev/null && pwd -P)/$(basename "$(readlink "$out1")")"
                resolved_out2="$(cd "$(dirname "$out2")" && cd "$(dirname "$(readlink "$out2")")" 2>/dev/null && pwd -P)/$(basename "$(readlink "$out2")")"

                resolved_in1="$(cd "$(dirname "$in1")" && pwd -P)/$(basename "$in1")"
                resolved_in2="$(cd "$(dirname "$in2")" && pwd -P)/$(basename "$in2")"

                [[ "$resolved_out1" == "$resolved_in1" && "$resolved_out2" == "$resolved_in2" ]] \
                    || die "Resume validation: unchanged Stage 03 symlink target mismatch for $lib"

            elif [[ "$action" == "downsample" ]]; then
                [[ ! -L "$out1" && ! -L "$out2" ]] \
                    || die "Resume validation: downsampled Stage 03 outputs unexpectedly symlinked for $lib"

            else
                die "Unknown normalisation action for $lib: $action"
            fi

            pass "resume: validated existing Stage 03 output for $lib"
            continue
        fi

        if [[ -e "$out1" || -e "$out2" ]]; then
            die "Resume validation: partial/incomplete Stage 03 output for $lib"
        fi

    else
        [[ ! -e "$out1" && ! -e "$out2" ]] \
            || die "Unexpected pre-existing normalised output for $lib"
    fi

    if [[ "$action" == "downsample" ]]; then
        "$RASUSA_BIN" reads \
            --coverage "$DEPTH_CAP_X" \
            --genome-size "$PLANNING_GENOME_SIZE_BP" \
            --seed "$DOWNSAMPLE_SEED" \
            -o "$out1" \
            -o "$out2" \
            "$in1" \
            "$in2"

    elif [[ "$action" == "unchanged" ]]; then
        ln -s "$(cd "$(dirname "$in1")" && pwd -P)/$(basename "$in1")" "$out1"
        ln -s "$(cd "$(dirname "$in2")" && pwd -P)/$(basename "$in2")" "$out2"

    else
        die "Unknown normalisation action for $lib: $action"
    fi

    gzip -t "$out1" "$out2"

    n1="$("$SEQKIT_BIN" stats -T "$out1" | awk -F'\t' 'NR==2{gsub(/,/,"",$4);print $4}')"
    n2="$("$SEQKIT_BIN" stats -T "$out2" | awk -F'\t' 'NR==2{gsub(/,/,"",$4);print $4}')"

    [[ "$n1" -eq "$n2" ]] \
        || die "Normalisation broke pairing for $lib: R1=$n1 R2=$n2"

done < <(
    python3 - "$PLAN" <<'PYPLAN'
import csv
import sys

with open(sys.argv[1], newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        print(
            r["library_id"],
            r["normalisation_action"],
            sep="\t",
        )
PYPLAN
)

cp "$PLAN" "${ANA}/04_normalised/depth_normalisation.tsv"

stamp "STAGE 03 COMPLETE"

echo
echo "=============================================================="
echo "PRODUCTION CANDIDATE CURRENTLY ENDS AFTER STAGE 03."
echo "No assembly/QC/typing/tree stages are present yet."
echo "=============================================================="

# ======================================================================
# STAGE 04: DE NOVO ASSEMBLY
# ======================================================================

stamp "STAGE 04: SPAdes isolate assembly for all libraries"

mkdir -p "${ANA}/05_assembly"

while IFS=$'\t' read -r lib
do
    stamp "Stage 04: $lib"

    d="${ANA}/05_assembly/${lib}"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        if [[ -s "${d}/contigs.fasta" ]]; then
            grep -q '^>' "${d}/contigs.fasta" \
                || die "Resume validation: Stage 04 contigs lack FASTA records for $lib"

            pass "resume: validated existing Stage 04 assembly for $lib"
            continue
        fi

        if [[ -e "$d" ]]; then
            die "Resume validation: partial/incomplete Stage 04 assembly for $lib"
        fi

    else
        [[ ! -e "$d" ]] \
            || die "Unexpected pre-existing assembly directory for $lib"
    fi

    mkdir -p "$d"

    "$SPADES_BIN" \
        --isolate \
        -1 "${ANA}/04_normalised/${lib}/${lib}_R1.fastq.gz" \
        -2 "${ANA}/04_normalised/${lib}/${lib}_R2.fastq.gz" \
        -o "${d}/spades" \
        -t "$SPADES_THREADS" \
        -m "$SPADES_MEMORY_GB"

    [[ -s "${d}/spades/contigs.fasta" ]] \
        || die "SPAdes produced no contigs for $lib"

    cp "${d}/spades/contigs.fasta" "${d}/contigs.fasta"

    [[ -s "${d}/contigs.fasta" ]] \
        || die "Assembly copy failed for $lib"

done < <(
    python3 -c 'import csv,sys; f=open(sys.argv[1],newline=""); [print(r["library_id"]) for r in csv.DictReader(f,delimiter="\t")]' "$MANIFEST"
)

stamp "STAGE 04 COMPLETE"

# ======================================================================
# STAGE 05: QUAST + CHECKM + RIGBY QC
# ======================================================================

stamp "STAGE 05: QUAST + CheckM v1 + Rigby QC"

QUAST_ROOT="${ANA}/06_quast"
CHECKM_ROOT="${ANA}/07_checkm"

mkdir -p "$QUAST_ROOT" "$CHECKM_ROOT/bins"

while IFS=$'\t' read -r lib
do
    stamp "Stage 05 QUAST: $lib"

    qd="${QUAST_ROOT}/${lib}"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        if [[ -s "${qd}/report.tsv" ]]; then
            python3 - "${qd}/report.tsv" "$lib" <<'PYRESUME05Q'
import csv
import sys

path, lib = sys.argv[1:]
metrics = {}

with open(path) as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) >= 2:
            metrics[row[0]] = row[1]

for required in ("# contigs", "N50", "Total length"):
    if required not in metrics:
        raise SystemExit(
            f"Resume validation: QUAST metric {required!r} missing for {lib}"
        )

print(f"PASS: resume QUAST metrics present for {lib}")
PYRESUME05Q

            pass "resume: validated existing QUAST output for $lib"
        else
            if [[ -e "$qd" ]]; then
                die "Resume validation: partial/incomplete QUAST output for $lib"
            fi
        fi

    else
        [[ ! -e "$qd" ]] \
            || die "Unexpected pre-existing QUAST directory for $lib"
    fi

    if [[ ! -s "${qd}/report.tsv" ]]; then
        mkdir -p "$qd"

    "$CONDA" run --no-capture-output -n "$QUAST_ENV" \
        quast.py \
            --min-contig 500 \
            --threads 2 \
            -o "$qd" \
            "${ANA}/05_assembly/${lib}/contigs.fasta"

        [[ -s "${qd}/report.tsv" ]] \
            || die "QUAST report missing for $lib"
    fi

    if [[ ! -e "${CHECKM_ROOT}/bins/${lib}.fna" ]]; then
        ln -s \
            "${ANA}/05_assembly/${lib}/contigs.fasta" \
            "${CHECKM_ROOT}/bins/${lib}.fna"
    fi

done < <(
    python3 -c 'import csv,sys; f=open(sys.argv[1],newline=""); [print(r["library_id"]) for r in csv.DictReader(f,delimiter="\t")]' "$MANIFEST"
)

if [[ "$RUN_MODE" == "--resume" && -s "${CHECKM_ROOT}/checkm_qa.tsv" ]]; then

    python3 - "${CHECKM_ROOT}/checkm_qa.tsv" "$MANIFEST" <<'PYRESUME05C'
import csv
import sys

checkm_path, manifest_path = sys.argv[1:]

with open(manifest_path, newline="") as f:
    expected = {
        r["library_id"]
        for r in csv.DictReader(f, delimiter="\t")
    }

with open(checkm_path, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

observed = []
for r in rows:
    key = (r.get("Bin Id") or r.get("Bin Id ") or "").strip()
    if key:
        observed.append(key)

if len(observed) != len(set(observed)):
    raise SystemExit(
        "Resume validation: duplicate CheckM Bin Id values"
    )

if set(observed) != expected:
    missing = sorted(expected - set(observed))
    extra = sorted(set(observed) - expected)
    raise SystemExit(
        f"Resume validation: CheckM library mismatch; "
        f"missing={missing}; extra={extra}"
    )

print(
    f"PASS: resume CheckM table contains all "
    f"{len(expected)} expected libraries"
)
PYRESUME05C

    pass "resume: validated existing combined CheckM output"

else
    if [[ "$RUN_MODE" == "--resume" && -e "${CHECKM_ROOT}/checkm_qa.tsv" ]]; then
        die "Resume validation: empty/incomplete CheckM QA table"
    fi

    rm -rf "${CHECKM_ROOT}/work"

    "$CONDA" run --no-capture-output -n "$CHECKM_ENV" \
        checkm lineage_wf \
            -x fna \
            -t "$CHECKM_THREADS" \
            --pplacer_threads 1 \
            --tab_table \
            -f "${CHECKM_ROOT}/checkm_qa.tsv" \
            "${CHECKM_ROOT}/bins" \
            "${CHECKM_ROOT}/work"

    [[ -s "${CHECKM_ROOT}/checkm_qa.tsv" ]] \
        || die "CheckM QA output missing"
fi

python3 "$HELPER" qc \
    --taxonomy "${ANA}/03_taxonomy/taxonomy_summary.tsv" \
    --quast "$QUAST_ROOT" \
    --checkm "${CHECKM_ROOT}/checkm_qa.tsv" \
    --max-contam "$MAX_CONTAMINATION" \
    --min-comp "$MIN_COMPLETENESS" \
    --max-contigs "$MAX_CONTIGS" \
    --min-n50 "$MIN_N50" \
    --min-len "$MIN_GENOME_LENGTH" \
    --max-len "$MAX_GENOME_LENGTH" \
    --output "${CHECKM_ROOT}/library_taxonomy_qc.tsv"

[[ "$(awk 'END{print NR-1}' "${CHECKM_ROOT}/library_taxonomy_qc.tsv")" -eq "$EXPECTED_LIBRARIES" ]] \
    || die "QC summary does not contain all $EXPECTED_LIBRARIES libraries"

stamp "STAGE 05 COMPLETE"

# ======================================================================
# STAGE 06: SALMONELLA-SPECIFIC TYPING
# ======================================================================

stamp "STAGE 06: SISTR + Achtman MLST + AMRFinderPlus"

TYPING_ROOT="${ANA}/09_typing"

mkdir -p "$TYPING_ROOT"

while IFS=$'\t' read -r lib
do
    stamp "Stage 06: $lib"

    td="${TYPING_ROOT}/${lib}"
    asm="${ANA}/05_assembly/${lib}/contigs.fasta"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        sistr_out="${td}/sistr.tsv"
        mlst_out="${td}/mlst.tsv"
        amr_out="${td}/amrfinder.tsv"

        if [[ -s "$sistr_out" && -s "$mlst_out" && -s "$amr_out" ]]; then
            python3 - "$sistr_out" "$mlst_out" "$amr_out" "$lib" <<'PYRESUME06'
import csv
import sys

sistr_path, mlst_path, amr_path, lib = sys.argv[1:]

with open(sistr_path, newline="") as f:
    sistr = list(csv.DictReader(f, delimiter="\t"))

if len(sistr) != 1:
    raise SystemExit(
        f"Resume validation: expected one SISTR row for {lib}"
    )

with open(mlst_path) as f:
    mlst_line = f.read().strip()

if not mlst_line:
    raise SystemExit(
        f"Resume validation: empty MLST output for {lib}"
    )

with open(amr_path, newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    try:
        header = next(reader)
    except StopIteration:
        raise SystemExit(
            f"Resume validation: empty AMRFinder output for {lib}"
        )

if "Element symbol" not in header:
    raise SystemExit(
        f"Resume validation: unexpected AMRFinder header for {lib}"
    )
PYRESUME06

            pass "resume: validated existing Stage 06 typing for $lib"
            continue
        fi

        if [[ -e "$td" ]]; then
            die "Resume validation: partial/incomplete Stage 06 typing for $lib"
        fi

    else
        [[ ! -e "$td" ]] \
            || die "Unexpected pre-existing typing directory for $lib"
    fi

    mkdir -p "$td"

    "$CONDA" run --no-capture-output -n "$SISTR_ENV" \
        sistr \
            --qc \
            --alleles-output "${td}/sistr_alleles.json" \
            --cgmlst-profiles "${td}/sistr_cgmlst.csv" \
            -f tab \
            -o "${td}/sistr" \
            "$asm"

    [[ -s "${td}/sistr.tab" ]] \
        || die "SISTR tab output missing for $lib"

    mv "${td}/sistr.tab" "${td}/sistr.tsv"

    PATH="$NCBI_BLAST_BIN_DIR:$(dirname "$ANY2FASTA_BIN"):$PATH" \
        "$MLST_BIN" "$asm" > "${td}/mlst.tsv"

    [[ -s "${td}/mlst.tsv" ]] \
        || die "MLST output missing for $lib"

    "$CONDA" run --no-capture-output -n "$AMRFINDER_ENV" \
        amrfinder \
            -n "$asm" \
            -d "$AMRFINDER_DB" \
            --threads 2 \
            --plus \
            -o "${td}/amrfinder.tsv"

    [[ -s "${td}/amrfinder.tsv" ]] \
        || die "AMRFinderPlus output missing for $lib"

done < <(
    python3 - "${CHECKM_ROOT}/library_taxonomy_qc.tsv" <<'PYTYPE'
import csv
import sys

with open(sys.argv[1], newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["salmonella_gate_pass"] == "yes":
            print(r["library_id"])
PYTYPE
)

python3 "$HELPER" typing \
    --qc "${CHECKM_ROOT}/library_taxonomy_qc.tsv" \
    --typing "$TYPING_ROOT" \
    --output "${TYPING_ROOT}/library_taxonomy_qc_typing.tsv"

[[ "$(awk 'END{print NR-1}' "${TYPING_ROOT}/library_taxonomy_qc_typing.tsv")" -eq "$EXPECTED_LIBRARIES" ]] \
    || die "Typing summary does not retain all $EXPECTED_LIBRARIES libraries"

stamp "STAGE 06 COMPLETE"

echo
echo "=============================================================="
echo "STAGE 06 BUILD COMPLETE; proceeding to Stage 07."
echo "=============================================================="

# ======================================================================
# STAGE 07: KNOWN TECHNICAL-REPLICATE SNP COMPARISON + RECONCILIATION
# ======================================================================

stamp "STAGE 07: MLW/CGRdeep technical-replicate reconciliation"

TECH_ROOT="${ANA}/10_technical_replicates"
RECON_ROOT="${ANA}/11_reconciliation"

mkdir -p "$TECH_ROOT" "$RECON_ROOT"

TECH_PLAN="${TECH_ROOT}/technical_replicate_snippy_plan.tsv"
TECH_LONG="${TECH_ROOT}/technical_replicate_pairwise_distances.tsv"
TECH_MATRIX="${TECH_ROOT}/technical_replicate_snp_matrix.tsv"

# ----------------------------------------------------------------------
# 07A. Build the 26-pair plan from the completed Stage 06 summary.
# The provisional reference is chosen only for mapping comparability;
# it is NOT yet the canonical genome.
# ----------------------------------------------------------------------

python3 - \
    "${TYPING_ROOT}/library_taxonomy_qc_typing.tsv" \
    "$TECH_PLAN" <<'PY07PLAN'
import csv
import sys
from collections import defaultdict

summary_path, output_path = sys.argv[1:]

with open(summary_path, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

by_bio = defaultdict(list)

for r in rows:
    by_bio[r["biological_isolate_id"]].append(r)

def score(r):
    return (
        float(r["checkm_contamination"]),
        -float(r["checkm_completeness"]),
        int(r["contig_count"]),
        -int(r["n50_bp"]),
        r["library_id"],
    )

plan = []

for bio, items in sorted(by_bio.items()):
    mlw = [r for r in items if r["sequencing_stream"] == "MLW"]
    cgr = [r for r in items if r["sequencing_stream"] == "CGRdeep"]

    if not mlw or not cgr:
        continue

    if len(mlw) != 1 or len(cgr) != 1:
        raise SystemExit(
            f"Unexpected within-stream multiplicity for cross-stream biological isolate {bio}: "
            f"MLW={len(mlw)} CGRdeep={len(cgr)}"
        )

    pair = mlw + cgr

    if all(r["rigby_qc_pass"] == "yes" for r in pair):
        reference = sorted(pair, key=score)[0]["library_id"]
        status = "READY_FOR_SNIPPY"
    else:
        reference = "NA"
        status = "NOT_BOTH_RIGBY_QC_PASS"

    plan.append({
        "biological_isolate_id": bio,
        "MLW_library": mlw[0]["library_id"],
        "CGRdeep_library": cgr[0]["library_id"],
        "MLW_rigby_qc_pass": mlw[0]["rigby_qc_pass"],
        "CGRdeep_rigby_qc_pass": cgr[0]["rigby_qc_pass"],
        "provisional_reference_library": reference,
        "status": status,
    })

if len(plan) != 26:
    raise SystemExit(
        f"Expected exactly 26 known MLW/CGRdeep biological pairs; found {len(plan)}"
    )

fields = [
    "biological_isolate_id",
    "MLW_library",
    "CGRdeep_library",
    "MLW_rigby_qc_pass",
    "CGRdeep_rigby_qc_pass",
    "provisional_reference_library",
    "status",
]

with open(output_path, "w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(plan)

print(
    f"technical replicate pairs={len(plan)} "
    f"ready_for_snippy={sum(r['status']=='READY_FOR_SNIPPY' for r in plan)} "
    f"not_both_qc_pass={sum(r['status']=='NOT_BOTH_RIGBY_QC_PASS' for r in plan)}"
)
PY07PLAN

# ----------------------------------------------------------------------
# Helper: translate any real project path into its /data Docker path.
# This deliberately resolves host symlinks first.
# ----------------------------------------------------------------------

to_container_path() {
    python3 - "$PROJECT" "$1" <<'PY07PATH'
import sys
from pathlib import Path

project = Path(sys.argv[1]).resolve()
path = Path(sys.argv[2]).resolve()

try:
    rel = path.relative_to(project)
except ValueError:
    raise SystemExit(
        f"Path is outside TiNTS project and cannot be mounted safely: {path}"
    )

print("/data/" + rel.as_posix())
PY07PATH
}

# ----------------------------------------------------------------------
# 07B. Run both technical representations against the SAME provisional
# reference and calculate their common-core SNP distance.
# ----------------------------------------------------------------------

while IFS=$'\t' read -r bio mlw cgr ref
do
    stamp "Stage 07 Snippy: $bio [$mlw vs $cgr; ref=$ref]"

    pairdir="${TECH_ROOT}/${bio}"
    pair_result="${pairdir}/pairwise_distance.tsv"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        if [[ -s "$pair_result" ]]; then
            python3 - "$pair_result" "$bio" "$mlw" "$cgr" <<'PY07RESUME'
import csv
import sys

path, bio, mlw, cgr = sys.argv[1:]

with open(path, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

if len(rows) != 1:
    raise SystemExit(
        f"Resume validation: expected one Stage 07 distance row for {bio}"
    )

r = rows[0]

if (
    r.get("biological_isolate_id") != bio
    or r.get("MLW_library") != mlw
    or r.get("CGRdeep_library") != cgr
):
    raise SystemExit(
        f"Resume validation: Stage 07 pair identity mismatch for {bio}"
    )

float(r["core_snp_distance"])
PY07RESUME

            pass "resume: validated existing Stage 07 comparison for $bio"
            continue
        fi

        if [[ -e "$pairdir" ]]; then
            die "Resume validation: partial/incomplete Stage 07 output for $bio"
        fi

    else
        [[ ! -e "$pairdir" ]] \
            || die "Unexpected pre-existing Stage 07 directory for $bio"
    fi

    mkdir -p "$pairdir"

    ref_host="${ANA}/05_assembly/${ref}/contigs.fasta"
    [[ -s "$ref_host" ]] || die "Missing provisional reference assembly: $ref"

    ref_container="$(to_container_path "$ref_host")"

    for lib in "$mlw" "$cgr"
    do
        r1_host="${ANA}/04_normalised/${lib}/${lib}_R1.fastq.gz"
        r2_host="${ANA}/04_normalised/${lib}/${lib}_R2.fastq.gz"

        [[ -s "$r1_host" && -s "$r2_host" ]] \
            || die "Missing normalised reads for technical replicate $lib"

        r1_container="$(to_container_path "$r1_host")"
        r2_container="$(to_container_path "$r2_host")"

        out_host="${pairdir}/${lib}"
        out_container="$(to_container_path "$out_host")"

        docker run --rm \
            --platform "$SNIPPY_PLATFORM" \
            -v "$PROJECT:/data" \
            "$SNIPPY_IMAGE" \
            snippy \
                --cpus 4 \
                --outdir "$out_container" \
                --ref "$ref_container" \
                --R1 "$r1_container" \
                --R2 "$r2_container"

        [[ -s "${out_host}/snps.txt" && -s "${out_host}/snps.vcf" ]] \
            || die "Snippy output incomplete for $lib"
    done

    pair_container="$(to_container_path "$pairdir")"

    docker run --rm \
        --platform "$SNIPPY_PLATFORM" \
        -v "$PROJECT:/data" \
        "$SNIPPY_IMAGE" \
        snippy-core \
            --ref "$ref_container" \
            --prefix "${pair_container}/core" \
            "${pair_container}/${mlw}" \
            "${pair_container}/${cgr}"

    [[ -s "${pairdir}/core.aln" ]] \
        || die "Snippy core alignment missing for $bio"

    "$CONDA" run --no-capture-output -n "$SNP_TOOLS_ENV" \
        snp-dists "${pairdir}/core.aln" \
        > "${pairdir}/core_snp_distances.tsv"

    python3 - \
        "${pairdir}/core_snp_distances.tsv" \
        "$bio" "$mlw" "$cgr" "$ref" "$pair_result" <<'PY07DIST'
import csv
import sys

matrix_path, bio, mlw, cgr, ref, output = sys.argv[1:]

with open(matrix_path) as f:
    rd = csv.reader(f, delimiter="\t")
    header = next(rd)
    matrix = {
        row[0]: {
            header[i]: row[i]
            for i in range(1, len(header))
        }
        for row in rd
    }

try:
    d1 = float(matrix[mlw][cgr])
    d2 = float(matrix[cgr][mlw])
except Exception as exc:
    raise SystemExit(
        f"Could not extract MLW/CGRdeep core SNP distance for {bio}: {exc}"
    )

if d1 != d2:
    raise SystemExit(
        f"Asymmetric SNP distance for {bio}: {d1} vs {d2}"
    )

row = {
    "biological_isolate_id": bio,
    "MLW_library": mlw,
    "CGRdeep_library": cgr,
    "provisional_reference_library": ref,
    "core_snp_distance": f"{d1:g}",
}

with open(output, "w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=row.keys(),
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerow(row)

print(
    f"{bio}: {mlw} vs {cgr} = {d1:g} common-core SNPs"
)
PY07DIST

done < <(
    python3 - "$TECH_PLAN" <<'PY07READY'
import csv
import sys

with open(sys.argv[1], newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["status"] == "READY_FOR_SNIPPY":
            print(
                r["biological_isolate_id"],
                r["MLW_library"],
                r["CGRdeep_library"],
                r["provisional_reference_library"],
                sep="\t",
            )
PY07READY
)

# ----------------------------------------------------------------------
# 07C. Aggregate pairwise distances and construct the sparse matrix
# required by the reconciliation helper.
# ----------------------------------------------------------------------

python3 - \
    "${TYPING_ROOT}/library_taxonomy_qc_typing.tsv" \
    "$TECH_PLAN" \
    "$TECH_ROOT" \
    "$TECH_LONG" \
    "$TECH_MATRIX" <<'PY07MATRIX'
import csv
import sys
from pathlib import Path

summary_path, plan_path, tech_root, long_out, matrix_out = sys.argv[1:]

with open(summary_path, newline="") as f:
    summary = list(csv.DictReader(f, delimiter="\t"))

with open(plan_path, newline="") as f:
    plan = list(csv.DictReader(f, delimiter="\t"))

qc_pass = sorted(
    r["library_id"]
    for r in summary
    if r["rigby_qc_pass"] == "yes"
)

distances = {}

long_rows = []

for r in plan:
    if r["status"] != "READY_FOR_SNIPPY":
        continue

    bio = r["biological_isolate_id"]
    p = Path(tech_root) / bio / "pairwise_distance.tsv"

    if not p.is_file():
        raise SystemExit(
            f"Missing completed technical-replicate distance for {bio}"
        )

    with p.open(newline="") as f:
        rr = list(csv.DictReader(f, delimiter="\t"))

    if len(rr) != 1:
        raise SystemExit(
            f"Expected one pairwise distance row for {bio}"
        )

    x = rr[0]
    mlw = x["MLW_library"]
    cgr = x["CGRdeep_library"]
    d = float(x["core_snp_distance"])

    distances[(mlw, cgr)] = d
    distances[(cgr, mlw)] = d
    long_rows.append(x)

long_fields = [
    "biological_isolate_id",
    "MLW_library",
    "CGRdeep_library",
    "provisional_reference_library",
    "core_snp_distance",
]

with open(long_out, "w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=long_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(long_rows)

with open(matrix_out, "w", newline="") as f:
    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow([""] + qc_pass)

    for a in qc_pass:
        row = [a]

        for b in qc_pass:
            if a == b:
                row.append("0")
            elif (a, b) in distances:
                row.append(f"{distances[(a,b)]:g}")
            else:
                row.append("NA")

        w.writerow(row)

print(
    f"QC-pass libraries in sparse matrix={len(qc_pass)}; "
    f"technical pairs with SNP distances={len(long_rows)}"
)
PY07MATRIX

# ----------------------------------------------------------------------
# 07D. Reconcile technical replicates and choose canonical genomes.
# ----------------------------------------------------------------------

python3 "$HELPER" reconcile \
    --summary "${TYPING_ROOT}/library_taxonomy_qc_typing.tsv" \
    --dist "$TECH_MATRIX" \
    --max-tech-snp "$TECHNICAL_REPLICATE_MAX_SNP" \
    --overrides "${PROJECT}/00_metadata/replicate_overrides.tsv" \
    --reconciliation "${RECON_ROOT}/technical_replicate_reconciliation.tsv" \
    --canonical "${RECON_ROOT}/canonical_isolates.tsv"

[[ "$(awk 'END{print NR-1}' "${RECON_ROOT}/technical_replicate_reconciliation.tsv")" -eq "$EXPECTED_BIOLOGICAL_ISOLATES" ]] \
    || die "Reconciliation table does not contain all $EXPECTED_BIOLOGICAL_ISOLATES biological isolates"

review_count="$(
    awk -F'\t' '
        NR==1 {
            for(i=1;i<=NF;i++) h[$i]=i
            next
        }
        $(h["status"])=="REVIEW_REQUIRED" {
            n++
        }
        END {
            print n+0
        }
    ' "${RECON_ROOT}/technical_replicate_reconciliation.tsv"
)"

if [[ "$review_count" -gt 0 ]]; then
    die "$review_count technical-replicate group(s) require manual review before canonical-genome analysis"
fi

[[ "$(awk 'END{print NR-1}' "${RECON_ROOT}/canonical_isolates.tsv")" -eq "$EXPECTED_BIOLOGICAL_ISOLATES" ]] \
    || die "Canonical table does not contain all $EXPECTED_BIOLOGICAL_ISOLATES biological isolates"

stamp "STAGE 07 COMPLETE"

echo
echo "=============================================================="
echo "PRODUCTION CANDIDATE CURRENTLY ENDS AFTER STAGE 07."
echo "No Prokka/Panaroo/final broad-tree stage is present yet."
echo "Technical-replicate SNP threshold: <= ${TECHNICAL_REPLICATE_MAX_SNP}"
echo "This threshold is NOT a transmission threshold."
echo "=============================================================="

# ======================================================================
# STAGE 08: CANONICAL-GENOME PROKKA ANNOTATION
# ======================================================================

stamp "STAGE 08: Prokka annotation of canonical Salmonella genomes only"

ANNOT_ROOT="${ANA}/08_annotation"
CANONICAL="${RECON_ROOT}/canonical_isolates.tsv"

[[ -s "$CANONICAL" ]] \
    || die "Canonical-isolate table missing after Stage 07"

mkdir -p "$ANNOT_ROOT"

while IFS=$'\t' read -r lib
do
    stamp "Stage 08: $lib"

    ad="${ANNOT_ROOT}/${lib}"
    gff="${ad}/${lib}.gff"

    if [[ "$RUN_MODE" == "--resume" ]]; then
        if [[ -s "$gff" ]]; then
            grep -q '^##gff-version' "$gff" \
                || die "Resume validation: malformed Prokka GFF for $lib"

            pass "resume: validated existing Stage 08 annotation for $lib"
            continue
        fi

        if [[ -e "$ad" ]]; then
            die "Resume validation: partial/incomplete Stage 08 annotation for $lib"
        fi

    else
        [[ ! -e "$ad" ]] \
            || die "Unexpected pre-existing Stage 08 annotation directory for $lib"
    fi

    mkdir -p "$ad"

    asm_host="${ANA}/05_assembly/${lib}/contigs.fasta"
    [[ -s "$asm_host" ]] || die "Canonical assembly missing for $lib"

    asm_container="$(to_container_path "$asm_host")"
    ad_container="$(to_container_path "$ad")"

    locustag="$(
        printf '%s' "$lib" |
        tr -cd '[:alnum:]' |
        cut -c1-18
    )"

    docker run --rm \
        --platform "$PROKKA_PLATFORM" \
        -v "$PROJECT:/data" \
        "$PROKKA_IMAGE" \
        prokka \
            --force \
            --cpus 2 \
            --outdir "$ad_container" \
            --prefix "$lib" \
            --locustag "$locustag" \
            "$asm_container"

    [[ -s "$gff" ]] \
        || die "Prokka GFF missing for $lib"

done < <(
    python3 - "$CANONICAL" <<'PY08'
import csv
import sys

with open(sys.argv[1], newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["rigby_qc_pass"] == "yes":
            print(r["library_id"])
PY08
)

stamp "STAGE 08 COMPLETE"

# ======================================================================
# STAGE 09: FINAL CANONICAL SALMONELLA BROAD TREE
# ======================================================================

stamp "STAGE 09: canonical Salmonella Panaroo + IQ-TREE 3 broad phylogeny"

TREE_ROOT="${ANA}/12_final_tree"
GFF_ROOT="${TREE_ROOT}/gffs"
PANAROO_ROOT="${TREE_ROOT}/panaroo"

if [[ "$RUN_MODE" == "--resume" && -s "${TREE_ROOT}/TiNTS_broad.treefile" ]]; then

    [[ -s "${PANAROO_ROOT}/core_gene_alignment.aln" ]] \
        || die "Resume validation: final core alignment missing"

    [[ -s "${TREE_ROOT}/core_snp_distances.tsv" ]] \
        || die "Resume validation: final SNP-distance matrix missing"

    pass "resume: validated existing Stage 09 broad-tree outputs"

else
    if [[ "$RUN_MODE" == "--resume" && -e "$TREE_ROOT" ]]; then
        die "Resume validation: partial/incomplete Stage 09 broad-tree output"
    fi

    if [[ "$RUN_MODE" != "--resume" && -e "$TREE_ROOT" ]]; then
        die "Unexpected pre-existing Stage 09 broad-tree directory"
    fi

    mkdir -p "$GFF_ROOT"

    # Copy, rather than symlink, canonical GFFs so the Docker container
    # never encounters host-absolute symlink targets.
    python3 - "$CANONICAL" "$ANNOT_ROOT" "$GFF_ROOT" <<'PY09GFF'
import csv
import shutil
import sys
from pathlib import Path

canonical_path, annot_root, gff_root = sys.argv[1:]

annot_root = Path(annot_root)
gff_root = Path(gff_root)

with open(canonical_path, newline="") as f:
    rows = [
        r for r in csv.DictReader(f, delimiter="\t")
        if r["rigby_qc_pass"] == "yes"
    ]

if len(rows) < 4:
    raise SystemExit(
        f"Need at least four canonical Salmonella genomes for broad tree; found {len(rows)}"
    )

for r in rows:
    lib = r["library_id"]
    src = annot_root / lib / f"{lib}.gff"
    dst = gff_root / f"{lib}.gff"

    if not src.is_file():
        raise SystemExit(f"Missing canonical Prokka GFF: {src}")

    shutil.copy2(src, dst)

print(f"canonical Salmonella GFFs copied={len(rows)}")
PY09GFF

    mapfile -t PANAROO_GFFS < <(
        find "$GFF_ROOT" \
            -maxdepth 1 \
            -type f \
            -name '*.gff' \
            -print |
        sort |
        while read -r p
        do
            to_container_path "$p"
        done
    )

    [[ "${#PANAROO_GFFS[@]}" -ge 4 ]] \
        || die "Insufficient canonical GFFs for Panaroo"

    panaroo_container="$(to_container_path "$PANAROO_ROOT")"

    docker run --rm \
        --platform "$PANAROO_PLATFORM" \
        -v "$PROJECT:/data" \
        "$PANAROO_IMAGE" \
        panaroo \
            -i "${PANAROO_GFFS[@]}" \
            -o "$panaroo_container" \
            --clean-mode moderate \
            --remove-invalid-genes \
            --refind-mode off \
            --alignment core \
            --aligner mafft \
            --core_threshold "$PANAROO_CORE_THRESHOLD" \
            --family_threshold "$PANAROO_IDENTITY" \
            -t "$PANAROO_THREADS"

    CORE="${PANAROO_ROOT}/core_gene_alignment.aln"

    [[ -s "$CORE" ]] \
        || die "Panaroo core alignment missing"

    "$CONDA" run --no-capture-output -n "$SNP_TOOLS_ENV" \
        snp-sites \
            -o "${TREE_ROOT}/core_snps.aln" \
            "$CORE"

    [[ -s "${TREE_ROOT}/core_snps.aln" ]] \
        || die "Core SNP alignment missing"

    "$CONDA" run --no-capture-output -n "$SNP_TOOLS_ENV" \
        snp-dists \
            "${TREE_ROOT}/core_snps.aln" \
        > "${TREE_ROOT}/core_snp_distances.tsv"

    [[ -s "${TREE_ROOT}/core_snp_distances.tsv" ]] \
        || die "Core SNP distance matrix missing"

    "$IQTREE3_BIN" \
        -s "$CORE" \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T "$IQTREE_THREADS" \
        --seed "$TREE_SEED" \
        --prefix "${TREE_ROOT}/TiNTS_broad"

    [[ -s "${TREE_ROOT}/TiNTS_broad.treefile" ]] \
        || die "IQ-TREE 3 broad tree missing"
fi

stamp "STAGE 09 COMPLETE"

echo
echo "=============================================================="
echo "PRODUCTION CANDIDATE CURRENTLY ENDS AFTER STAGE 09."
echo "Broad tree uses one canonical Salmonella genome per biological isolate."
echo "No species-wide Gubbins analysis is performed."
echo "=============================================================="

# ======================================================================
# STAGE 10: TREE RENDERING, CANDIDATES, LEDGER, RESULTS EXPORT
# ======================================================================

stamp "STAGE 10: final interpretation outputs and RESULTS area"

CANDIDATE_ROOT="${ANA}/13_finescale_candidates"
mkdir -p "$CANDIDATE_ROOT"

# ----------------------------------------------------------------------
# 10A. Render broad tree and prepare tree metadata.
# Household is retained in metadata but is not used as a broad-tree track.
# ----------------------------------------------------------------------

"$CONDA" run --no-capture-output -n "$PLOT_ENV" \
    python "$HELPER" plot-tree \
        --tree "${TREE_ROOT}/TiNTS_broad.treefile" \
        --canonical "$CANONICAL" \
        --metadata "$ISOLATE_METADATA" \
        --tree-metadata "${TREE_ROOT}/tree_metadata.tsv" \
        --svg "${TREE_ROOT}/TiNTS_Salmonella_broad_tree.svg" \
        --pdf "${TREE_ROOT}/TiNTS_Salmonella_broad_tree.pdf"

[[ -s "${TREE_ROOT}/tree_metadata.tsv" ]] \
    || die "Final tree metadata missing"

[[ -s "${TREE_ROOT}/TiNTS_Salmonella_broad_tree.svg" ]] \
    || die "Final broad-tree SVG missing"

[[ -s "${TREE_ROOT}/TiNTS_Salmonella_broad_tree.pdf" ]] \
    || die "Final broad-tree PDF missing"

# ----------------------------------------------------------------------
# 10B. Candidate groups for later lineage/reference-specific analysis.
# This is candidate identification only; no direct transmission inference.
# ----------------------------------------------------------------------

python3 "$HELPER" candidates \
    --canonical "$CANONICAL" \
    --metadata "$ISOLATE_METADATA" \
    --output "${CANDIDATE_ROOT}/finescale_candidate_groups.tsv"

[[ -s "${CANDIDATE_ROOT}/finescale_candidate_groups.tsv" ]] \
    || die "Fine-scale candidate-group table missing"

# ----------------------------------------------------------------------
# 10C. Explicit exclusion/disposition ledger.
# No sequencing library should silently disappear.
# ----------------------------------------------------------------------

python3 "$HELPER" ledger \
    --summary "${TYPING_ROOT}/library_taxonomy_qc_typing.tsv" \
    --reconciliation "${RECON_ROOT}/technical_replicate_reconciliation.tsv" \
    --output "${RECON_ROOT}/exclusion_ledger.tsv"

[[ -s "${RECON_ROOT}/exclusion_ledger.tsv" ]] \
    || die "Exclusion ledger missing"

[[ "$(awk 'END{print NR-1}' "${RECON_ROOT}/exclusion_ledger.tsv")" -eq "$EXPECTED_LIBRARIES" ]] \
    || die "Exclusion ledger does not contain all $EXPECTED_LIBRARIES libraries"

# ----------------------------------------------------------------------
# 10D. RESULTS collision guard.
# RESULTS is written only after the whole analytical workflow succeeds.
# ----------------------------------------------------------------------

if output_has_content "$RESULTS_DIR"; then
    die "RESULTS directory unexpectedly contains files before final export"
fi

mkdir -p "$RESULTS_DIR"
RESULTS_ENV="${RESULTS_DIR}/environment"
mkdir -p "$RESULTS_ENV"

# ----------------------------------------------------------------------
# 10E. Primary tabular outputs.
# ----------------------------------------------------------------------

cp "$MANIFEST" \
    "${RESULTS_DIR}/TiNTS_library_manifest.tsv"

cp "${ANA}/02_host_depleted/host_depletion_summary.tsv" \
    "${RESULTS_DIR}/TiNTS_host_depletion_summary.tsv"

cp "${ANA}/03_taxonomy/taxonomy_summary.tsv" \
    "${RESULTS_DIR}/TiNTS_taxonomy.tsv"

cp "${ANA}/04_normalised/depth_normalisation.tsv" \
    "${RESULTS_DIR}/TiNTS_depth_normalisation.tsv"

cp "${TYPING_ROOT}/library_taxonomy_qc_typing.tsv" \
    "${RESULTS_DIR}/TiNTS_library_taxonomy_qc_typing.tsv"

cp "${RECON_ROOT}/exclusion_ledger.tsv" \
    "${RESULTS_DIR}/TiNTS_exclusion_ledger.tsv"

cp "${RECON_ROOT}/technical_replicate_reconciliation.tsv" \
    "${RESULTS_DIR}/TiNTS_technical_replicates.tsv"

cp "$CANONICAL" \
    "${RESULTS_DIR}/TiNTS_canonical_isolates.tsv"

cp "${CANDIDATE_ROOT}/finescale_candidate_groups.tsv" \
    "${RESULTS_DIR}/TiNTS_finescale_candidate_groups.tsv"

cp "$TECH_PLAN" \
    "${RESULTS_DIR}/TiNTS_technical_replicate_snippy_plan.tsv"

cp "$TECH_LONG" \
    "${RESULTS_DIR}/TiNTS_technical_replicate_pairwise_distances.tsv"

# ----------------------------------------------------------------------
# 10F. Broad-tree outputs.
# ----------------------------------------------------------------------

cp "${TREE_ROOT}/TiNTS_broad.treefile" \
    "${RESULTS_DIR}/TiNTS_Salmonella_broad_tree.nwk"

cp "${TREE_ROOT}/TiNTS_Salmonella_broad_tree.svg" \
    "${RESULTS_DIR}/TiNTS_Salmonella_broad_tree.svg"

cp "${TREE_ROOT}/TiNTS_Salmonella_broad_tree.pdf" \
    "${RESULTS_DIR}/TiNTS_Salmonella_broad_tree.pdf"

cp "${TREE_ROOT}/tree_metadata.tsv" \
    "${RESULTS_DIR}/TiNTS_Salmonella_tree_metadata.tsv"

cp "${PANAROO_ROOT}/core_gene_alignment.aln" \
    "${RESULTS_DIR}/TiNTS_Salmonella_core_alignment.fasta"

cp "${TREE_ROOT}/core_snp_distances.tsv" \
    "${RESULTS_DIR}/TiNTS_Salmonella_core_snp_distances.tsv"

# ----------------------------------------------------------------------
# 10G. Provenance.
# ----------------------------------------------------------------------

cp "$INPUT_CHECKSUMS" \
    "${RESULTS_DIR}/input_checksums.sha256"

cp "$SOFTWARE_VERSIONS" \
    "${RESULTS_DIR}/software_versions.tsv"

cp "$DATABASE_VERSIONS" \
    "${RESULTS_DIR}/database_versions.tsv"

cp "$CHECKM_ENV_EXPORT" "$RESULTS_ENV/"
cp "$QUAST_ENV_EXPORT" "$RESULTS_ENV/"
cp "$SISTR_ENV_EXPORT" "$RESULTS_ENV/"
cp "$AMRFINDER_ENV_EXPORT" "$RESULTS_ENV/"
cp "$SNP_TOOLS_ENV_EXPORT" "$RESULTS_ENV/"
cp "$GUBBINS_ENV_EXPORT" "$RESULTS_ENV/"
cp "$PLOT_ENV_EXPORT" "$RESULTS_ENV/"
cp "$MLST_ENV_FILE" "$RESULTS_ENV/"
cp "$BLAST_ENV_FILE" "$RESULTS_ENV/"

cp "${SCRIPT_DIR}/tints_config.env" \
    "${RESULTS_ENV}/tints_config.env"

cp "${BASH_SOURCE[0]}" \
    "${RESULTS_ENV}/run_tints_genomics.sh"

cp "$HELPER" \
    "${RESULTS_ENV}/tints_helpers.py"

{
    echo -e "container\timage\tplatform\tdigest"
    echo -e "Prokka\t${PROKKA_IMAGE}\t${PROKKA_PLATFORM}\t${PROKKA_IMAGE_DIGEST}"
    echo -e "Panaroo\t${PANAROO_IMAGE}\t${PANAROO_PLATFORM}\t${PANAROO_IMAGE_DIGEST}"
    echo -e "Snippy\t${SNIPPY_IMAGE}\t${SNIPPY_PLATFORM}\t${SNIPPY_IMAGE_DIGEST}"
} > "${RESULTS_ENV}/container_images.tsv"

shasum -a 256 \
    "${RESULTS_ENV}/tints_config.env" \
    "${RESULTS_ENV}/run_tints_genomics.sh" \
    "${RESULTS_ENV}/tints_helpers.py" \
    > "${RESULTS_ENV}/workflow_checksums.sha256"

# ----------------------------------------------------------------------
# 10H. Results README.
# ----------------------------------------------------------------------

{
    echo "TiNTS GENOMICS RESULTS"
    echo "Generated: $(date -Iseconds)"
    echo
    echo "SCOPE"
    echo "- Descriptive bacterial isolate genomics supporting the TiNTS epidemiological manuscript."
    echo "- 111 sequencing libraries represent 85 biological isolates."
    echo "- MLW and CGRdeep representations remain separate until explicit technical-replicate reconciliation."
    echo
    echo "HOST-SEQUENCE SCREENING"
    echo "- Every library, including environmental isolates, was screened against ${HUMAN_REFERENCE_NAME} (${HUMAN_REFERENCE_ACCESSION})."
    echo "- fastp preprocessing preceded Bowtie2 ${HOST_FILTER_PRESET} host screening."
    echo "- A read pair was retained only when BOTH mates were unmapped to the human reference."
    echo "- Kraken2 subsequently provided an independent residual-human taxid ${HUMAN_TAXID} audit."
    echo "- Computational screening does not prove literal absence of every human-derived fragment."
    echo
    echo "TAXONOMY AND DEPTH"
    echo "- Kraken2 taxonomy was performed on the full host-depleted library before any depth downsampling."
    echo "- Kraken2 production runs used the 2026-02-26 Standard-16 database in RAM without --memory-mapping."
    echo "- Salmonella clade abundance >=${SALMONELLA_READ_PCT_THRESHOLD}% was required for Salmonella-specific analysis."
    echo "- Non-Salmonella and mixed/ambiguous libraries remain explicitly represented in the output tables."
    echo "- Libraries above ${DEPTH_CAP_X}x planning depth were deterministically capped with Rasusa using seed ${DOWNSAMPLE_SEED}; lower-depth libraries were unchanged."
    echo
    echo "ASSEMBLY AND QC"
    echo "- SPAdes --isolate was used for all libraries."
    echo "- QUAST and CheckM v1 provided general assembly/QC metrics."
    echo "- Rigby Salmonella QC thresholds: contamination <=${MAX_CONTAMINATION}%; completeness >=${MIN_COMPLETENESS}%; contigs <=${MAX_CONTIGS}; N50 >=${MIN_N50} bp; genome length ${MIN_GENOME_LENGTH}-${MAX_GENOME_LENGTH} bp inclusive."
    echo
    echo "SALMONELLA TYPING"
    echo "- Salmonella-gate-pass libraries underwent SISTR, Achtman MLST, and AMRFinderPlus."
    echo
    echo "TECHNICAL REPLICATES"
    echo "- Known MLW/CGRdeep representations were compared only after independent taxonomy, assembly, QC and typing."
    echo "- Both representations were mapped to the same provisional QC-selected reference for technical-resequencing comparison."
    echo "- <=${TECHNICAL_REPLICATE_MAX_SNP} common-core SNPs is used only as a stringent technical-replicate consistency criterion."
    echo "- This SNP criterion is NOT a transmission threshold."
    echo "- Canonical genome choice is QC-based, not sequencing-depth-based."
    echo
    echo "BROAD SALMONELLA TREE"
    echo "- One canonical Salmonella genome per biological isolate was annotated with Prokka after reconciliation."
    echo "- Panaroo core-presence threshold: ${PANAROO_CORE_THRESHOLD}; family identity threshold: ${PANAROO_IDENTITY}."
    echo "- The broad phylogeny used native IQ-TREE 3 with ModelFinder, 1000 ultrafast bootstraps and 1000 SH-aLRT replicates."
    echo "- Gubbins was NOT applied species-wide across the diverse broad Salmonella tree."
    echo
    echo "INTERPRETATION"
    echo "- Genomic similarity does not prove direct person-to-person or environment-to-person transmission."
    echo "- Fine-scale genomic analysis requires lineage/reference-specific follow-up, common callable-genome comparison, temporal/household context and relevant background diversity."
    echo
    echo "PRIMARY FILES"
    find "$RESULTS_DIR" -maxdepth 1 -type f -exec basename {} \; | sort
} > "${RESULTS_DIR}/README_RESULTS.txt"

# ----------------------------------------------------------------------
# 10I. Final completeness audit.
# ----------------------------------------------------------------------

for required in \
    README_RESULTS.txt \
    TiNTS_library_manifest.tsv \
    TiNTS_host_depletion_summary.tsv \
    TiNTS_taxonomy.tsv \
    TiNTS_depth_normalisation.tsv \
    TiNTS_library_taxonomy_qc_typing.tsv \
    TiNTS_exclusion_ledger.tsv \
    TiNTS_technical_replicates.tsv \
    TiNTS_canonical_isolates.tsv \
    TiNTS_finescale_candidate_groups.tsv \
    TiNTS_Salmonella_broad_tree.nwk \
    TiNTS_Salmonella_broad_tree.svg \
    TiNTS_Salmonella_broad_tree.pdf \
    TiNTS_Salmonella_tree_metadata.tsv \
    TiNTS_Salmonella_core_alignment.fasta \
    TiNTS_Salmonella_core_snp_distances.tsv \
    input_checksums.sha256 \
    software_versions.tsv \
    database_versions.tsv
do
    [[ -s "${RESULTS_DIR}/${required}" ]] \
        || die "Required final RESULTS file missing/empty: $required"
done

stamp "STAGE 10 COMPLETE"
stamp "TiNTS GENOMICS PRODUCTION WORKFLOW COMPLETE"

echo
echo "=============================================================="
echo "TiNTS genomics workflow completed successfully."
echo "Results: ${RESULTS_DIR}"
echo
echo "Genomic similarity must not be interpreted as proof of direct transmission."
echo "=============================================================="
