# TiNTS bacterial genomics pipeline

This directory contains the bacterial-genomics workflows used for the TiNTS study.

## Contents

- `scripts/run_tints_genomics.sh` — primary production pipeline
- `scripts/run_tints_salmonella_rescue.sh` — mixed-library Salmonella rescue workflow
- `scripts/tints_helpers.py` — helper functions used by the production pipeline
- `scripts/tints_config.env` — analysis parameters, software locations, QC thresholds, and database configuration
- `environments/` — frozen Conda environment specifications
- `provenance/software_versions.tsv` — software versions used for the analysis
- `provenance/database_versions.tsv` — reference/database versions and checksums
- `examples/example_library_manifest.tsv` — example input-manifest structure

## Primary workflow

Reads are processed with fastp, depleted of human sequence against T2T-CHM13v2.0 using Bowtie2, taxonomically classified with Kraken2 Standard-16, depth-normalised with Rasusa, assembled with SPAdes, and assessed with QUAST and CheckM.

Salmonella genomes are typed using SISTR and Achtman MLST and characterised with AMRFinderPlus. Downstream steps include technical-replicate reconciliation, annotation, core-genome analysis, and phylogenetic inference.

The main entry point is:

```bash
bash scripts/run_tints_genomics.sh
```

## Mixed-library Salmonella rescue

The rescue workflow extracts reads assigned by Kraken2 to *Salmonella* taxid 590 and descendants, requires adequate Salmonella-specific sequencing depth, normalises target depth, and reassembles those reads.

Rescued assemblies are required to pass the same final genome-quality criteria used in the primary analysis; rescue-entry criteria are triage criteria and do not represent relaxed final QC.

Entry point:

```bash
bash scripts/run_tints_salmonella_rescue.sh
```

## Configuration

The production analysis used the parameters recorded in `scripts/tints_config.env`.

By default the scripts expect the working genomics project at:

```text
~/TiNTS_genomics
```

A different project root can be supplied without editing the frozen configuration:

```bash
export TINTS_PROJECT=/path/to/TiNTS_genomics
```

If Conda is installed somewhere other than `~/miniforge3/bin/conda`, set:

```bash
export TINTS_CONDA=/path/to/conda
```

The pipeline expects the metadata, reads, reference databases, and analysis directories described in the configuration file. The example manifest illustrates the required library-level metadata structure.

## Software and reference databases

Software versions used for the manuscript analysis are recorded in `provenance/software_versions.tsv`.

Database/reference snapshots and checksums are recorded in `provenance/database_versions.tsv`, including:

- Kraken2 Standard-16
- CheckM
- AMRFinderPlus
- T2T-CHM13v2.0 human host-depletion reference

The Conda YAML files under `environments/` have machine-specific installation prefixes removed so that they can be recreated on another system.

## Genome-quality criteria

The final assembly QC thresholds used by the production workflow are recorded in `scripts/tints_config.env`, including completeness, contamination, contig number, N50, and genome-length thresholds.

## Data availability

Raw sequencing reads and large intermediate files are not stored in this Git repository. Sequencing data and assembled genomes are deposited separately in the European Nucleotide Archive as described in the manuscript.

## Reproducibility

`SHA256SUMS.txt` records checksums for the files in this directory. The production scripts retain fixed seeds where deterministic behaviour is required and record the relevant software/database versions and analysis parameters.
