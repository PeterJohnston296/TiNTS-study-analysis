# TiNTS study analysis

This repository contains the R/Quarto workflows used for the epidemiological analyses supporting the TiNTS manuscript.

## R analysis workflow

Run the documents in order:

1. `01_variable_preparation.qmd`
2. `02_elastic_net_selection.qmd`
3. `03_cox_models.qmd`
4. `04_episode_burden_and_poisson.qmd`
5. `05_ct_threshold_sensitivity.qmd`

The scripts expect the corresponding analysis inputs under:

- `data_raw/`
- `data_clean/`
- `rds/`

Generated outputs are written beneath `outputs/`.

The public workflows are streamlined from the development analysis notebooks. Exploratory code, superseded analyses, and manuscript-layout code have been omitted, while the analysis steps required to fit the reported models and generate the principal numerical outputs are retained.

## Bayesian models

The hidden Markov models are supplied separately with their Stan code and R/Quarto calling workflows.

The full Bayesian fits are computationally intensive, particularly Model B, and may require multi-day runtime. The repository therefore documents the model specification, data-preparation workflow, priors, and sampler settings, together with the fitted posterior summaries used for the manuscript where appropriate.

## Data availability

Participant-level analytic data are not included in this public repository.

A pseudonymised analytic dataset sufficient to run the R workflows can be made available for confidential peer review and, subject to the study's ethical approvals and data-sharing arrangements, through controlled access for research use.

The public repository does not contain the participant, household, or sample-level analytic datasets.

## Reproducibility

The R workflows are designed to be run from the repository root with the corresponding analytic inputs placed in `data_raw/`, `data_clean/`, and `rds/`.

Package requirements are specified within the individual Quarto documents. Analysis outputs are written to the `outputs/` directory.
