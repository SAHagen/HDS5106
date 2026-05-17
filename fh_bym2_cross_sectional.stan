data {
  int<lower=1> N;
  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=N> obs_idx;
  
  vector[N_obs] y;
  vector<lower=0>[N_obs] psi;
  
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  real<lower=0> scaling_factor;
  
  int<lower=0> K;
  matrix[N, K] X;
}

parameters {
  real mu;
  vector[N] v;
  vector[N] u_raw;
  real<lower=0> sigma_b;
  real<lower=0, upper=1> phi;
  vector[K] beta;
}

transformed parameters {
  vector[N] u = u_raw - mean(u_raw);
  vector[N] b = sigma_b * (sqrt(1 - phi) * v + sqrt(phi / scaling_factor) * u);
  
  vector[N] theta;
  for (i in 1:N) {
    theta[i] = mu + b[i];
  }
  if (K > 0) {
    theta += X * beta;
  }
}

model {
  mu ~ normal(-1.85, 1);
  sigma_b ~ exponential(4.6);
  phi ~ beta(3, 2);
  beta ~ normal(0, 1);
  v ~ normal(0, 1);
  
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);
  
  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  vector<lower=0, upper=1>[N] prevalence;
  for (i in 1:N) {
    prevalence[i] = inv_logit(theta[i]);
  }
  real p_national = inv_logit(mu);
}