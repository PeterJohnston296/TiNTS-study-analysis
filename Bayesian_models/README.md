# Bayesian hidden Markov models

This directory contains the Bayesian models used in the TiNTS manuscript.

## Model A

`01_model_A.qmd` fits the individual two-state susceptible–infected hidden Markov model to the daily PCR observation grid.

Stan program:

`stan/model_A_individual_hmm.stan`

Required analytic input:

`data_clean/HMM_ModelA_input_dailygrid.csv`

## Model B

`02_model_B.qmd` fits the household-coupled hidden Markov model. Acquisition probability is modelled as a function of preceding household latent infection pressure, baseline household environmental contamination, baseline mobility, and a household random effect.

Stan program:

`stan/model_B_household_chmm.stan`

Required controlled-access analytic inputs include:

- `rds/model_B_prepared_input.rds`
- `data_raw/CRF_B_MCA.csv`

The household-coupled model has a large joint latent state space and is computationally intensive. On the development hardware, the full Model B fit required approximately three days. The scripts are therefore supplied principally to document and permit exact refitting of the model; routine review does not require rerunning the full fit.

Sampler settings, priors, random seeds, and threading configuration are recorded directly in the Quarto documents.
