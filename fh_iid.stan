
data {
  int<lower=1> N;
  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N> obs_idx;
  vector[N_obs] y;
  vector<lower=0>[N_obs] psi;
  int<lower=1> K;
  matrix[N, K] X;
}
parameters {
  real mu;
  vector[N] v;
  real<lower=0> sigma_b;
  vector[K] beta;
}
transformed parameters {
  vector[N] theta;
  theta = mu + sigma_b * v + X * beta;
}
model {
  mu ~ normal(-1.85, 1.0);
  sigma_b ~ exponential(4.6);
  beta ~ normal(0, 1);
  v ~ std_normal();
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

