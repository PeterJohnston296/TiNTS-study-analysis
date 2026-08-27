# TiNTS study analysis

This repository contains the analysis code supporting the TiNTS manuscript:

**Frequent asymptomatic carriage and household transmission of non-typhoidal Salmonella in urban Malawi: a prospective longitudinal cohort and modelling study**

The repository contains three main components.

## 1. Epidemiological analyses

The top-level Quarto documents reproduce the principal R-based epidemiological analyses:

1. `01_variable_preparation.qmd`
2. `02_elastic_net_selection.qmd`
3. `03_cox_models.qmd`
4. `04_episode_burden_and_poisson.qmd`
5. `05_ct_threshold_sensitivity.qmd`

These cover preparation of the epidemiological analysis dataset, elastic-net variable selection, first-event and recurrent-event Cox models, episode burden and recurrence analyses, and Ct-threshold sensitivity analyses.

The workflows expect the corresponding analytic inputs under:

- `data_raw/`
- `data_clean/`
- `rds/`

Generated outputs are written beneath `outputs/`.

The public workflows are streamlined from the development analysis notebooks. Exploratory code, superseded analyses, and manuscript-layout code have been omitted while retaining the analysis steps required to fit the reported models and generate the principal numerical outputs.

## 2. Bayesian transmission models

The `bayesian_models/` directory contains the hidden Markov models used to estimate latent acquisition, clearance, and household dependence.

It includes:

- **Model A:** individual susceptible–infected hidden Markov model;
- **Model B:** household-coupled hidden Markov model;
- the corresponding Stan model files;
- R/Quarto workflows documenting data preparation, priors, sampler settings, model fitting, and posterior summaries.

Model B is computationally intensive because household latent states are modelled jointly and may require multi-day runtime depending on hardware.

See `bayesian_models/README.md` for details.

## 3. Genomics pipeline

The `genomics_pipeline/` directory contains the bacterial genomics workflow used for the TiNTS sequencing analyses.

It includes:

- the primary production Bash pipeline;
- the mixed-library *Salmonella* rescue workflow;
- helper code and frozen analysis configuration;
- Conda environment specifications;
- software and database provenance;
- an example library manifest;
- SHA256 checksums.

The workflow includes read QC, human-read depletion, taxonomic classification, assembly, genome QC, *Salmonella* typing, AMR characterisation, annotation, core-genome analysis, and phylogenetic inference.

See `genomics_pipeline/README.md` for details.

## Data availability

Participant-level analytic data are not included in this public repository.

A pseudonymised analytic dataset sufficient to run the R workflows can be made available for confidential peer review and, subject to the study's ethical approvals and data-sharing agreements, through controlled access for research use.

Raw sequencing reads and assembled genomes are deposited separately in the European Nucleotide Archive as described in the manuscript.

Study standard operating procedures, case record forms, participant information documents, and consent and assent materials are deposited separately on Zenodo.

## Reproducibility

Package requirements and analysis settings are documented within the individual workflows.

The Bayesian and genomics directories contain additional model, software, database, and configuration information required to reproduce those components of the study.

This repository represents the analysis-code version corresponding to the submitted manuscript.
