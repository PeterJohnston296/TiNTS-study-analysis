// -----------------------------------------------------------------------------
// Prior-sensitivity version of Model A
//
// This model keeps the same likelihood and latent-state structure as Model A,
// but uses alternative priors to test robustness.
//
// Key points:
//   1. Two-state SIS hidden Markov model:
//        S = susceptible / not latently infected
//        I = latently infected
//
//   2. Ct-aware two-stage emission model:
//        sampled day with no Ct      -> informative non-amplification
//        sampled day with Ct value   -> amplification + state-specific Ct density
//        unsampled day               -> no emission contribution
//
//   3. The likelihood is evaluated using the forward algorithm.
//
//   4. Generated quantities return smoothed posterior infection probabilities:
//        prob_I[n] = P(Z_n = I | all observed data for that participant)
//      not merely filtered probabilities.
// -----------------------------------------------------------------------------

functions {
  vector emission_log_probs(
    int has_obs,
    int has_ct,
    real ct,
    real phi_S,
    real phi_I,
    real mu_S,
    real mu_I,
    real sigma_S,
    real sigma_I
  ) {
    vector[2] log_em;

    // State order: S = 1, I = 2

    if (has_obs == 0) {
      // No PCR test was performed: emission contribution is 1.
      log_em[1] = 0;
      log_em[2] = 0;

    } else {
      if (has_ct == 1) {
        // PCR amplification occurred and a Ct value was observed.
        log_em[1] = log(phi_S) + normal_lpdf(ct | mu_S, sigma_S);
        log_em[2] = log(phi_I) + normal_lpdf(ct | mu_I, sigma_I);

      } else {
        // PCR test was performed but no amplification occurred.
        log_em[1] = log1m(phi_S);
        log_em[2] = log1m(phi_I);
      }
    }

    return log_em;
  }
}

data {
  // ---------------------------------------------------------
  // INDEXING
  // ---------------------------------------------------------
  int<lower=1> N;                                // total person-days
  int<lower=1> N_subj;                           // number of participants

  array[N_subj] int<lower=1, upper=N> start_idx; // first row for each participant
  array[N_subj] int<lower=1, upper=N> end_idx;   // final row for each participant

  // ---------------------------------------------------------
  // OBSERVATION MODEL INPUTS
  // ---------------------------------------------------------
  // has_obs[n] = 0: no sample/PCR test performed
  // has_obs[n] = 1 and has_ct[n] = 0: PCR performed, no amplification
  // has_obs[n] = 1 and has_ct[n] = 1: PCR performed, Ct value observed

  array[N] int<lower=0, upper=1> has_obs;
  array[N] int<lower=0, upper=1> has_ct;

  vector[N] Ct;              // only used when has_ct[n] == 1
  real<lower=0> Ct_max;      // retained for compatibility; not used in this model
}

parameters {
  // ---------------------------------------------------------
  // LATENT STATE DYNAMICS
  // ---------------------------------------------------------
  real alpha_SI;    // logit daily acquisition probability
  real alpha_IS;    // logit daily clearance probability
  real alpha_pi;    // logit initial latent infection probability

  // ---------------------------------------------------------
  // PCR AMPLIFICATION GATE
  // ---------------------------------------------------------
  real alpha_phi_I; // logit P(amplification | infected)
  real alpha_phi_S; // logit P(amplification | susceptible)

  // ---------------------------------------------------------
  // Ct DISTRIBUTIONS GIVEN AMPLIFICATION
  // ---------------------------------------------------------
  // Constraint enforces that Ct values in the susceptible state are higher
  // than those in the infected state.

  real<lower=10, upper=43.5> mu_I;
  real<lower=mu_I + 1.5, upper=45> mu_S;

  real<lower=1, upper=6> sigma_I;
  real<lower=1, upper=6> sigma_S;
}

transformed parameters {
  // ---------------------------------------------------------
  // PROBABILITY SCALE PARAMETERS
  // ---------------------------------------------------------
  real<lower=0, upper=1> p_SI  = inv_logit(alpha_SI);
  real<lower=0, upper=1> p_IS  = inv_logit(alpha_IS);
  real<lower=0, upper=1> pi_I  = inv_logit(alpha_pi);

  real<lower=0, upper=1> phi_I = inv_logit(alpha_phi_I);
  real<lower=0, upper=1> phi_S = inv_logit(alpha_phi_S);

  // ---------------------------------------------------------
  // TRANSITION MATRIX ON LOG SCALE
  // ---------------------------------------------------------
  // State order: S = 1, I = 2
  // Rows are previous state; columns are current state.

  matrix[2,2] log_T;

  log_T[1,1] = log1m(p_SI);  // S -> S
  log_T[1,2] = log(p_SI);    // S -> I
  log_T[2,1] = log(p_IS);    // I -> S
  log_T[2,2] = log1m(p_IS);  // I -> I
}

model {
  // ---------------------------------------------------------
  // PRIOR-SENSITIVITY PRIORS
  // ---------------------------------------------------------
  // These are deliberately different from the primary Model A priors:
  // lower acquisition, longer infection duration, fewer false-positive
  // amplifications, and less certain detection during latent infection.

  alpha_SI    ~ normal(logit(0.01), 0.7);
  alpha_IS    ~ normal(logit(0.17), 0.5);
  alpha_pi    ~ normal(logit(0.063), 0.4);

  alpha_phi_I ~ normal(logit(0.75), 0.35);
  alpha_phi_S ~ normal(logit(0.005), 0.8);

  mu_I        ~ normal(28, 2);
  mu_S        ~ normal(40, 2);

  sigma_I     ~ normal(4, 0.5);
  sigma_S     ~ normal(3, 0.5);

  // ---------------------------------------------------------
  // LIKELIHOOD VIA FORWARD ALGORITHM
  // ---------------------------------------------------------
  // The first row for each participant receives the initial state probability
  // directly. There is no transition before the first modelled time point.

  for (s in 1:N_subj) {
    int start = start_idx[s];
    int end   = end_idx[s];

    vector[2] log_alpha_prev;
    vector[2] log_alpha_curr;
    vector[2] log_em;

    // First time point: initial distribution plus first emission.
    log_em = emission_log_probs(
      has_obs[start],
      has_ct[start],
      Ct[start],
      phi_S,
      phi_I,
      mu_S,
      mu_I,
      sigma_S,
      sigma_I
    );

    log_alpha_prev[1] = log1m(pi_I) + log_em[1];
    log_alpha_prev[2] = log(pi_I)  + log_em[2];

    // Subsequent time points: transition then emission.
    if (end > start) {
      for (n in (start + 1):end) {
        log_em = emission_log_probs(
          has_obs[n],
          has_ct[n],
          Ct[n],
          phi_S,
          phi_I,
          mu_S,
          mu_I,
          sigma_S,
          sigma_I
        );

        log_alpha_curr[1] =
          log_em[1] +
          log_sum_exp(
            log_alpha_prev[1] + log_T[1,1],
            log_alpha_prev[2] + log_T[2,1]
          );

        log_alpha_curr[2] =
          log_em[2] +
          log_sum_exp(
            log_alpha_prev[1] + log_T[1,2],
            log_alpha_prev[2] + log_T[2,2]
          );

        log_alpha_prev = log_alpha_curr;
      }
    }

    // Marginal likelihood contribution for this participant.
    target += log_sum_exp(log_alpha_prev);
  }
}

generated quantities {
  // ---------------------------------------------------------
  // SMOOTHED POSTERIOR INFECTION PROBABILITIES
  // ---------------------------------------------------------
  // prob_I[n] is now:
  //   P(Z_n = I | all observed data for that participant)
  //
  // This is forward-backward smoothing, not filtering.

  vector[N] prob_I;
  vector[N_subj] log_lik_subj;

  for (s in 1:N_subj) {
    int start = start_idx[s];
    int end   = end_idx[s];
    int L     = end - start + 1;

    array[L] vector[2] log_alpha;
    array[L] vector[2] log_beta;

    vector[2] log_em;
    vector[2] log_smooth;

    // -------------------------------------------------------
    // Forward pass
    // -------------------------------------------------------

    log_em = emission_log_probs(
      has_obs[start],
      has_ct[start],
      Ct[start],
      phi_S,
      phi_I,
      mu_S,
      mu_I,
      sigma_S,
      sigma_I
    );

    log_alpha[1][1] = log1m(pi_I) + log_em[1];
    log_alpha[1][2] = log(pi_I)  + log_em[2];

    if (L > 1) {
      for (ell in 2:L) {
        int n = start + ell - 1;

        log_em = emission_log_probs(
          has_obs[n],
          has_ct[n],
          Ct[n],
          phi_S,
          phi_I,
          mu_S,
          mu_I,
          sigma_S,
          sigma_I
        );

        log_alpha[ell][1] =
          log_em[1] +
          log_sum_exp(
            log_alpha[ell - 1][1] + log_T[1,1],
            log_alpha[ell - 1][2] + log_T[2,1]
          );

        log_alpha[ell][2] =
          log_em[2] +
          log_sum_exp(
            log_alpha[ell - 1][1] + log_T[1,2],
            log_alpha[ell - 1][2] + log_T[2,2]
          );
      }
    }

    log_lik_subj[s] = log_sum_exp(log_alpha[L]);

    // -------------------------------------------------------
    // Backward pass
    // -------------------------------------------------------
    // log_beta[ell][z] is the probability of future observations
    // after time ell, conditional on state z at time ell.

    log_beta[L][1] = 0;
    log_beta[L][2] = 0;

    if (L > 1) {
      for (rev_ell in 1:(L - 1)) {
        int ell = L - rev_ell;
        int n_next = start + ell;

        vector[2] log_em_next;

        log_em_next = emission_log_probs(
          has_obs[n_next],
          has_ct[n_next],
          Ct[n_next],
          phi_S,
          phi_I,
          mu_S,
          mu_I,
          sigma_S,
          sigma_I
        );

        log_beta[ell][1] =
          log_sum_exp(
            log_T[1,1] + log_em_next[1] + log_beta[ell + 1][1],
            log_T[1,2] + log_em_next[2] + log_beta[ell + 1][2]
          );

        log_beta[ell][2] =
          log_sum_exp(
            log_T[2,1] + log_em_next[1] + log_beta[ell + 1][1],
            log_T[2,2] + log_em_next[2] + log_beta[ell + 1][2]
          );
      }
    }

    // -------------------------------------------------------
    // Smoothed posterior probabilities
    // -------------------------------------------------------

    for (ell in 1:L) {
      int n = start + ell - 1;

      log_smooth[1] = log_alpha[ell][1] + log_beta[ell][1];
      log_smooth[2] = log_alpha[ell][2] + log_beta[ell][2];

      prob_I[n] = exp(log_smooth[2] - log_sum_exp(log_smooth));
    }
  }
}
