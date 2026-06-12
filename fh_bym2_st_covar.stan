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
  //real<lower=0> sigma_delta;                        // temporal SD
  
  vector[N] gamma_raw;  // raw interaction
  real<lower=0> sigma_gamma;
  vector[K] beta;                                   // covariate coefficients
}

transformed parameters {
   // BYM2 district effect
  vector[N] u = u_raw - mean(u_raw);
  vector[N] b = sigma_b * (sqrt(1 - phi) * v + sqrt(phi / scaling_factor) * u);

  // Sum-to-zero temporal effect
  vector[T] delta;
  delta[1] =  alpha;
  delta[2] = -alpha;

  // Space-time interaction, sum-to-zero across districts within round 1
  // Antisymmetric across time: γ_{i,2} = -γ_{i,1}
  vector[N] gamma_star = gamma_raw - mean(gamma_raw);

  // Linear predictor
  // Linear predictor
  vector[NT] theta;
  for (n in 1:NT) {
    real interaction;
    if (time[n] == 1) {
      interaction =  sigma_gamma * gamma_star[district[n]];
    } else {
      interaction = -sigma_gamma * gamma_star[district[n]];
    }
    theta[n] = mu + b[district[n]] + delta[time[n]] + interaction;
  }
  // Covariate contribution
  if (K > 0) {
    theta += X * beta;
  }
}

model {
  // Priors
  mu          ~ normal(-2, 1);
  alpha       ~ normal(0, 2.5);     // direct weakly informative on temporal effect
  sigma_b     ~ exponential(4.6);
  phi         ~ beta(3, 2);
  sigma_gamma ~ exponential(2);     // SD of trajectory deviations
  beta        ~ normal(0, 1);

  // Non-centred parameterisations
  v         ~ std_normal();
  gamma_raw ~ std_normal();

  // ICAR prior + soft sum-to-zero
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);

  // Likelihood
  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  // District-round prevalence
  vector<lower=0, upper=1>[NT] prevalence;
  for (n in 1:NT) prevalence[n] = inv_logit(theta[n]);

  // National change on logit scale (signed: negative = decline)
  real delta_change = -2 * alpha;

  // National-level (intercept-only) round prevalences
  real p_round1 = inv_logit(mu + alpha);
  real p_round2 = inv_logit(mu - alpha);

  // District-level trajectory deviation on logit scale
  // (negative = faster decline than national; positive = slower / rising)
  vector[N] district_trend_dev;
  for (i in 1:N) {
    district_trend_dev[i] = -2 * sigma_gamma * gamma_star[i];
  }
}

