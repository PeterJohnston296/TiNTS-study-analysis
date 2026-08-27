# TiNTS manuscript R reproducibility code

This repository contains the R/Quarto workflows supporting the epidemiological analyses reported in the TiNTS manuscript.

## Analysis order

1. `01_variable_preparation.qmd`
2. `02_elastic_net_selection.qmd`
3. `03_cox_models.qmd`
4. `04_episode_burden_and_poisson.qmd`
5. `05_ct_threshold_sensitivity.qmd`

The first two documents form a sequential workflow. The remaining documents reconstruct the relevant analysis inputs from the de-identified/pseudonymised analytic data package.

## Data

Participant-level source data are **not stored in this GitHub repository**. The analysis expects the following directory structure at the repository/project root:

```text
data_raw/
  CRF_A.csv
  CRF_B_MCA.csv
  household_distance_from_river.csv
data_clean/
  wealth_index_by_household.csv
  bayesian_pcr_for_epi_analysis_28052025.csv
  pcr_stool_id_lookup.csv
  ct_threshold_sensitivity/
rds/
  episodes_df.rds
  followup_df.rds
```

The exact analytic input package should be distributed only under the study's approved data-sharing arrangements. Study materials and SOPs are archived separately at Zenodo (doi:10.5281/zenodo.15690901).

## Reproducibility

Run from the repository root with Quarto/R, for example:

```bash
quarto render 01_variable_preparation.qmd
quarto render 02_elastic_net_selection.qmd
quarto render 03_cox_models.qmd
quarto render 04_episode_burden_and_poisson.qmd
quarto render 05_ct_threshold_sensitivity.qmd
```

Outputs are written beneath `outputs/`.
