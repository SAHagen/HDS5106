
data {
  int<lower=1> N;
  int<lower=1> C;
  array[N] int<lower=1, upper=C> country;

  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N> obs_idx;
  vector[N_obs] y;
  vector<lower=0>[N_obs] psi;

  int<lower=1> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  real<lower=0> scaling_factor;

  int<lower=1> K;
  matrix[N, K] X;
}
parameters {
  real mu;
  vector[C] alpha_raw;
  vector[N] v;
  vector[N] u_raw;
  real<lower=0> sigma_b;
  real<lower=0, upper=1> phi;
  vector[K] beta_bar;
  vector<lower=0>[K] tau;
  matrix[C, K] beta_z;
}
transformed parameters {
  vector[C] alpha = alpha_raw - mean(alpha_raw);
  vector[N] u = u_raw - mean(u_raw);
  vector[N] b = sigma_b * (sqrt(1 - phi) * v
                           + sqrt(phi / scaling_factor) * u);
  matrix[C, K] beta;
  for (c in 1:C)
    for (k in 1:K)
      beta[c, k] = beta_bar[k] + tau[k] * beta_z[c, k];

  vector[N] theta;
  for (i in 1:N) {
    int ci = country[i];
    theta[i] = mu + alpha[ci] + b[i] + dot_product(X[i], beta[ci]');
  }
}
model {
  mu        ~ normal(-1.85, 1.0);
  alpha_raw ~ normal(0, 1.0);
  sigma_b   ~ exponential(4.6);
  phi       ~ beta(3, 2);

  v         ~ std_normal();
  target   += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);

  beta_bar  ~ normal(0, 1);
  tau       ~ exponential(2);
  to_vector(beta_z) ~ std_normal();

  y ~ normal(theta[obs_idx], psi);
}
generated quantities {
  vector[N_obs] log_lik;
  for (i in 1:N_obs)
    log_lik[i] = normal_lpdf(y[i] | theta[obs_idx[i]], psi[i]);

  vector[N] prevalence;
  for (i in 1:N)
    prevalence[i] = inv_logit(theta[i]);

  vector[C] p_national;
  for (c in 1:C) {
    real sum_p = 0;
    int  n_c   = 0;
    for (i in 1:N) {
      if (country[i] == c) {
        sum_p += prevalence[i];
        n_c   += 1;
      }
    }
    p_national[c] = sum_p / n_c;
  }
}

