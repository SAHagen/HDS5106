
data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1> NT;
  int<lower=1> C;
  int<lower=1> K;

  int<lower=1> N_obs;
  array[N_obs] int<lower=1, upper=NT> obs_idx;
  vector[N_obs] y;
  vector<lower=0>[N_obs] psi;

  array[NT] int<lower=1, upper=N> district;
  array[NT] int<lower=1, upper=T> time;
  vector[NT] year_c;
  array[N] int<lower=1, upper=C> country;

  matrix[NT, K] X;

  int<lower=1> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  vector<lower=0>[C] scaling_factor_country;
}

parameters {
  real mu;
  vector[C - 1] alpha_raw;

  vector[C] delta;
  real delta_bar;
  real<lower=0> sigma_delta;

  vector[N] v;
  vector[N] u_raw;
  vector<lower=0>[C] sigma_b;
  vector<lower=0, upper=1>[C] phi;

  matrix[C, K] beta_raw;
  vector[K] beta_bar;
  vector<lower=0>[K] tau;

  vector[NT] gamma_raw;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  vector[C] alpha;
  alpha[1] = 0;
  alpha[2:C] = alpha_raw;

  matrix[C, K] beta;
  for (k in 1:K)
    for (c in 1:C)
      beta[c, k] = beta_bar[k] + tau[k] * beta_raw[c, k];

  vector[N] u = u_raw - mean(u_raw);

  vector[N] b;
  for (i in 1:N) {
    int c = country[i];
    b[i] = sigma_b[c] * ( sqrt(1 - phi[c]) * v[i]
                          + sqrt(phi[c] / scaling_factor_country[c]) * u[i] );
  }

  vector[NT] gamma_star = sigma_gamma * gamma_raw;

  vector[NT] theta;
  for (n in 1:NT) {
    int d = district[n];
    int c = country[d];

    real xb = 0;
    for (k in 1:K)
      xb += X[n, k] * beta[c, k];

    theta[n] = mu + alpha[c] + b[d]
               + delta[c] * year_c[n]
               + gamma_star[n]
               + xb;
  }
}

model {
  mu ~ normal(-1.85, 1.5);
  alpha_raw ~ normal(0, 1.5);

  delta_bar ~ normal(0, 0.5);
  sigma_delta ~ exponential(5);
  delta ~ normal(delta_bar, sigma_delta);

  sigma_b ~ exponential(4.6);
  phi ~ beta(3, 2);
  v ~ std_normal();
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);

  beta_bar ~ normal(0, 1);
  tau ~ exponential(2);
  to_vector(beta_raw) ~ std_normal();

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

  matrix[C, T] p_national;
  for (c in 1:C) {
    for (t in 1:T) {
      real total = 0;
      int count = 0;
      for (n in 1:NT) {
        if (country[district[n]] == c && time[n] == t) {
          total += prevalence[n];
          count += 1;
        }
      }
      p_national[c, t] = total / count;
    }
  }
}

