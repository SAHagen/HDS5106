
data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1> NT;
  int<lower=1> K;

  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=NT> obs_idx;
  vector[N_obs] y;
  vector<lower=0>[N_obs] psi;

  array[NT] int<lower=1, upper=N> district;
  array[NT] int<lower=1, upper=T> time;
  vector[NT] year_c;

  matrix[NT, K] X;

  int<lower=1> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  real<lower=0> scaling_factor;
}

parameters {
  real mu;
  real delta;

  vector[N] v;
  vector[N] u_raw;
  real<lower=0> sigma_b;
  real<lower=0, upper=1> phi;

  vector[K] beta;

  vector[NT] gamma_raw;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  vector[N] u = u_raw - mean(u_raw);
  vector[N] b = sigma_b * (sqrt(1 - phi) * v
                           + sqrt(phi / scaling_factor) * u);

  vector[NT] gamma_star = sigma_gamma * gamma_raw;

  vector[NT] theta;
  for (n in 1:NT) {
    int d = district[n];
    real xb = 0;
    for (k in 1:K)
      xb += X[n, k] * beta[k];
    theta[n] = mu + b[d] + delta * year_c[n] + gamma_star[n] + xb;
  }
}

model {
  mu ~ normal(-1.85, 1.5);
  delta ~ normal(0, 0.5);

  sigma_b ~ exponential(4.6);
  phi ~ beta(3, 2);
  v ~ std_normal();
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);

  beta ~ normal(0, 1);

  sigma_gamma ~ exponential(2);
  gamma_raw ~ std_normal();

  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  vector[N_obs] log_lik;
  for (i in 1:N_obs)
    log_lik[i] = normal_lpdf(y[i] | theta[obs_idx[i]], psi[i]);

  vector[NT] prevalence;
  for (n in 1:NT)
    prevalence[n] = inv_logit(theta[n]);

  // National prevalence at each round (population-weighted mean would be
  // ideal; here we use the unweighted district mean for consistency with
  // Model 2 of the workbook).
  vector[T] p_national;
  for (t in 1:T) {
    real total = 0;
    int  count = 0;
    for (n in 1:NT) {
      if (time[n] == t) { total += prevalence[n]; count += 1; }
    }
    p_national[t] = total / count;
  }

  // Round 1 and Round 2 calendar years (centred) and the implied
  // round-to-round shift on the logit scale.
  real yr1_c = -999;
  real yr2_c = -999;
  for (n in 1:NT) {
    if (time[n] == 1 && yr1_c == -999) yr1_c = year_c[n];
    if (time[n] == 2 && yr2_c == -999) yr2_c = year_c[n];
  }
  real delta_R2_R1     = delta * (yr2_c - yr1_c);
  real odds_ratio_R2_R1 = exp(delta_R2_R1);
}

