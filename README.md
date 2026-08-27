# TiNTS study analysis

This repository contains the R/Quarto workflows used for the epidemiological analyses supporting the TiNTS manuscript.

## R analysis workflow

Run the documents in order:

1. `01_variable_preparation.qmd`
2. `02_elastic_net_selection.qmd`
3. `03_cox_models.qmd`
4. `04_episode_burden_and_poisson.qmd`
5. `05_ct_threshold_sensitivity.qmd`

The scripts expect the analysis inputs to be available under:

- `data_raw/`
- `data_clean/`
- `rds/`

Outputs are written beneath `outputs/`.

The public workflows are streamlined from the development notebooks: exploratory code, superseded analyses, and manuscript-layout code have been omitted, while the analysis steps required to fit the reported models and generate the principal numerical outputs are retained.

## Bayesian models

The hidden Markov models are supplied separately with their Stan code and R/Quarto calling workflows. The full Bayesian fits are computationally intensive, particularly Model B, so the repository also preserves the fitted posterior summaries used for the manuscript where appropriate.

## Data

Participant-level analytic data are not committed to the public repository. Access and sharing are subject to the study's approved data-governance arrangements.
