functions {

  real household_loglik(
    int h,
    int T,
    array[] int S_h,
    array[] int M_h,
    array[] int X_h,
    array[,] int z_state,
    array[,,] int has_obs,
    array[,,] int has_ct,
    array[,,] real Ct,
    array[,] real D_star,
    real pi_I,
    real p_IS,
    real phi_I,
    real phi_S,
    real alpha_SI,
    real beta_env,
    real beta_hh,
    real beta_comm,
    vector u,
    real mu_I,
    real mu_S,
    real sigma_I,
    real sigma_S
  ) {
    int M = M_h[h];
    int S = S_h[h];

    vector[S] log_alpha_prev;
    vector[S] log_alpha_curr;

    // Precompute # infected in each joint state
    array[S] int Icount;
    for (s in 1:S) {
      int c = 0;
      for (m in 1:M) c += z_state[s, m];
      Icount[s] = c;
    }

    // Initial distribution factorises across members (prior over Z_{t=1})
    for (s in 1:S) {
      real lp = 0;
      for (m in 1:M) {
        int z = z_state[s, m];
        lp += (z == 1) ? log(pi_I) : log1m(pi_I);
      }
      log_alpha_prev[s] = lp;
    }

    // ---- t = 1: multiply by emissions only (NO transition yet) ----
    {
      vector[S] log_em;
      for (s in 1:S) {
        real le = 0;
        for (m in 1:M) {
          int z = z_state[s, m];
          if (has_obs[h, m, 1] == 1) {
            if (has_ct[h, m, 1] == 1) {
              if (z == 1)
                le += log(phi_I) + normal_lpdf(Ct[h, m, 1] | mu_I, sigma_I);
              else
                le += log(phi_S) + normal_lpdf(Ct[h, m, 1] | mu_S, sigma_S);
            } else {
              le += (z == 1) ? log1m(phi_I) : log1m(phi_S);
            }
          }
        }
        log_em[s] = le;
      }
      log_alpha_prev = log_alpha_prev + log_em;
    }

    // ---- t = 2..T: transition then emission ----
    for (t in 2:T) {
      vector[S] log_em;

      // Emissions at time t
      for (s in 1:S) {
        real le = 0;
        for (m in 1:M) {
          int z = z_state[s, m];
          if (has_obs[h, m, t] == 1) {
            if (has_ct[h, m, t] == 1) {
              if (z == 1)
                le += log(phi_I) + normal_lpdf(Ct[h, m, t] | mu_I, sigma_I);
              else
                le += log(phi_S) + normal_lpdf(Ct[h, m, t] | mu_S, sigma_S);
            } else {
              le += (z == 1) ? log1m(phi_I) : log1m(phi_S);
            }
          }
        }
        log_em[s] = le;
      }

      // Transition + forward update
      for (s2 in 1:S) {
        real acc = negative_infinity();

        for (s1 in 1:S) {
          int k_prev = Icount[s1];
          real lt = 0;

          for (m in 1:M) {
            int z1 = z_state[s1, m];
            int z2 = z_state[s2, m];

            real p_SI_m = inv_logit(alpha_SI
                                    + beta_env * X_h[h]
                                    + u[h]
                                    + beta_hh * k_prev
                                    + beta_comm * D_star[h, m]);

            if (z1 == 1)
              lt += (z2 == 0) ? log(p_IS) : log1m(p_IS);
            else
              lt += (z2 == 1) ? log(p_SI_m) : log1m(p_SI_m);
          }

          acc = log_sum_exp(acc, log_alpha_prev[s1] + lt);
        }

        log_alpha_curr[s2] = log_em[s2] + acc;
      }

      log_alpha_prev = log_alpha_curr;
    }

    return log_sum_exp(log_alpha_prev);
  }

  // ---- reduce_sum partial_sum (SLICE-AWARE indexing) ----
  real partial_sum(
    array[] int hh_slice, int start, int end,
    int T,
    array[] int S_h,
    array[] int M_h,
    array[] int X_h,
    array[,] int z_state,
    array[,,] int has_obs,
    array[,,] int has_ct,
    array[,,] real Ct,
    array[,] real D_star,
    real pi_I,
    real p_IS,
    real phi_I,
    real phi_S,
    real alpha_SI,
    real beta_env,
    real beta_hh,
    real beta_comm,
    vector u,
    real mu_I,
    real mu_S,
    real sigma_I,
    real sigma_S
  ) {
    real lp = 0;
    int N = end - start + 1;

    for (i in 1:N) {
      int h = hh_slice[i];
      lp += household_loglik(
        h, T, S_h, M_h, X_h, z_state, has_obs, has_ct, Ct, D_star,
        pi_I, p_IS, phi_I, phi_S,
        alpha_SI, beta_env, beta_hh, beta_comm, u,
        mu_I, mu_S, sigma_I, sigma_S
      );
    }
    return lp;
  }
}

data {
  int<lower=1> H;
  int<lower=1> M_max;
  int<lower=1> T;

  int<lower=1> S_max;
  array[H] int<lower=1> S_h;
  array[S_max, M_max] int<lower=0,upper=1> z_state;

  array[H] int<lower=1, upper=M_max> M_h;
  array[H] int<lower=0, upper=1> X_h;

  array[H, M_max, T] int<lower=0, upper=1> has_obs;
  array[H, M_max, T] int<lower=0, upper=1> has_ct;
  array[H, M_max, T] real Ct;

  array[H, M_max] real D_star;

  int<lower=1> grainsize;
}

parameters {
  real alpha_SI;
  real beta_env;
  real beta_hh;
  real beta_comm;

  real alpha_IS;
  real alpha_pi;

  real<lower=0> sigma_u;
  vector[H] z_u;

  real alpha_phi_I;
  real alpha_phi_S;

  ordered[2] mu_raw;

  real<lower=1, upper=6> sigma_I;
  real<lower=1, upper=6> sigma_S;
}

transformed parameters {
  vector[H] u = sigma_u * z_u;

  real mu_I = 10 + 35 * inv_logit(mu_raw[1]);
  real mu_S = 10 + 35 * inv_logit(mu_raw[2]);
}

model {
  alpha_SI ~ normal(logit(0.04), 0.8);
  alpha_IS ~ normal(logit(0.30), 0.7);
  alpha_pi ~ normal(logit(0.063), 0.4);

  beta_env  ~ normal(0, 0.7);
  beta_hh   ~ normal(0, 0.3);
  beta_comm ~ normal(0, 0.5);

  z_u     ~ normal(0, 1);
  sigma_u ~ normal(0, 0.4);

  alpha_phi_I ~ normal(logit(0.85), 0.25);
  alpha_phi_S ~ normal(logit(0.03), 0.8);

  mu_I ~ normal(28, 2);
  mu_S ~ normal(38.5, 2);

  sigma_I ~ normal(4, 0.7);
  sigma_S ~ normal(3, 0.7);

  {
    real pi_I  = inv_logit(alpha_pi);
    real p_IS  = inv_logit(alpha_IS);
    real phi_I = inv_logit(alpha_phi_I);
    real phi_S = inv_logit(alpha_phi_S);

    array[H] int hh_ids;
    for (h in 1:H) hh_ids[h] = h;

    target += reduce_sum(
      partial_sum,
      hh_ids,
      grainsize,
      T, S_h, M_h, X_h, z_state, has_obs, has_ct, Ct, D_star,
      pi_I, p_IS, phi_I, phi_S,
      alpha_SI, beta_env, beta_hh, beta_comm, u,
      mu_I, mu_S, sigma_I, sigma_S
    );
  }
}

generated quantities {
  vector[H] log_lik_household;

  {
    real pi_I  = inv_logit(alpha_pi);
    real p_IS  = inv_logit(alpha_IS);
    real phi_I = inv_logit(alpha_phi_I);
    real phi_S = inv_logit(alpha_phi_S);

    vector[H] u_local = sigma_u * z_u;

    for (h in 1:H) {
      log_lik_household[h] = household_loglik(
        h, T, S_h, M_h, X_h, z_state, has_obs, has_ct, Ct, D_star,
        pi_I, p_IS, phi_I, phi_S,
        alpha_SI, beta_env, beta_hh, beta_comm, u_local,
        mu_I, mu_S, sigma_I, sigma_S
      );
    }
  }
}
