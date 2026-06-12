
data {
  int<lower=1> N;
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
  vector[N] u_raw;
  real<lower=0> sigma_b;
  vector[K] beta;
}
transformed parameters {
  vector[N] theta;
  vector[N] u = u_raw - mean(u_raw);
  theta = mu + sigma_b * u / sqrt(scaling_factor) + X * beta;
}
model {
  mu ~ normal(-1.85, 1.0);
  sigma_b ~ exponential(4.6);
  beta ~ normal(0, 1);
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);
  y ~ normal(theta[obs_idx], psi);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N] prevalence;
  for (i in 1:N_obs)
    log_lik[i] = normal_lpdf(y[i] | theta[obs_idx[i]], psi[i]);
  for (i in 1:N)
    prevalence[i] = inv_logit(theta[i]);
}

