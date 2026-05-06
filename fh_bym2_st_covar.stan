//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// =============================================================================
// Fay-Herriot BYM2 Space-Time Model WITH COVARIATES
// Single country, two survey rounds
// Extends fh_bym2_st_null.stan with district-level static covariates
// =============================================================================

data {
  int<lower=1> N;                                   // number of districts
  int<lower=1> T;                                   // number of time points (= 2)
  int<lower=1> NT;                                  // total rows (N × T)

  int<lower=1> N_obs;                               // number of observed rows
  array[N_obs] int<lower=1, upper=NT> obs_idx;      // indices of observed rows

  array[NT] int<lower=1, upper=N> district;         // district index per row
  array[NT] int<lower=1, upper=T> time;             // time index per row

  vector[N_obs] y;                                  // observed logit direct estimates
  vector<lower=0>[N_obs] psi;                       // sampling SEs

  int<lower=0> N_edges;                             // spatial edges
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;

  real<lower=0> scaling_factor;                     // BYM2 scaling

  // ===== COVARIATES =====
  int<lower=0> K;                                   // number of covariates
  matrix[NT, K] X;                                  // covariate matrix
}

parameters {
  real mu;                                          // national intercept
  real alpha;                                       // temporal magnitude

  vector[N] v;                                      // IID district effects
  vector[N] u_raw;                                  // unscaled ICAR effects

  real<lower=0> sigma_b;                            // district SD
  real<lower=0, upper=1> phi;                       // BYM2 mixing
  real<lower=0> sigma_delta;                        // temporal SD

  vector[K] beta;                                   // covariate coefficients
}

transformed parameters {
  // Sum-to-zero ICAR
  vector[N] u = u_raw - mean(u_raw);

  // BYM2 district effect
  vector[N] b = sigma_b * (sqrt(1 - phi) * v + sqrt(phi / scaling_factor) * u);

  // Sum-to-zero temporal effect
  vector[T] delta;
  delta[1] =  alpha;
  delta[2] = -alpha;

  // Linear predictor
  vector[NT] theta;
  for (n in 1:NT) {
    theta[n] = mu + b[district[n]] + delta[time[n]];
  }

  // Covariate contribution
  if (K > 0) {
    theta += X * beta;
  }
}

model {
  // Priors
  mu          ~ normal(-2, 1);
  sigma_b     ~ exponential(4.6);
  phi         ~ beta(3, 2);
  sigma_delta ~ exponential(2);
  alpha       ~ normal(0, sigma_delta);

  v ~ std_normal();

  // ICAR prior
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);

  // Covariate coefficients (weakly informative for standardised X)
  beta ~ normal(0, 1);

  // Likelihood
  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  vector<lower=0, upper=1>[NT] prevalence;
  for (n in 1:NT) {
    prevalence[n] = inv_logit(theta[n]);
  }

  real delta_change = 2 * alpha;
  real p_round1 = inv_logit(mu + alpha);
  real p_round2 = inv_logit(mu - alpha);
}

