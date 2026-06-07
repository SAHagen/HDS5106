// =============================================================================
// fh_bym2_cross_sectional.stan
//
// Bayesian Fay-Herriot Small Area Estimation model with BYM2 spatial random
// effect, for a single survey round.
//
// References:
//   - Wakefield, Okonek & Pedersen (2020): Fay-Herriot framework
//   - Riebler et al. (2016): BYM2 reparameterisation with scaling factor
//   - Besag, York & Mollié (1991): Original BYM specification
//
// Inputs:
//   y       — logit-transformed direct prevalence estimates (observed districts)
//   psi     — design-based standard errors on the logit scale (delta method)
//   node1, node2 — undirected edges of the district adjacency graph
//   scaling_factor — Riebler et al. (2016) BYM2 scaling factor
//   X       — N x K covariate matrix (standardised)
// =============================================================================

data {
  int<lower=1> N;                                  // number of districts
  int<lower=1> N_obs;                              // number of observed districts
  array[N_obs] int<lower=1, upper=N> obs_idx;      // indices of observed districts
  vector[N_obs] y;                                 // logit direct estimates
  vector<lower=0>[N_obs] psi;                      // logit-scale SEs

  // Spatial adjacency graph
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  real<lower=0> scaling_factor;

  // Covariates
  int<lower=0> K;
  matrix[N, K] X;
}

parameters {
  real mu;                                         // national logit intercept
  vector[N] v;                                     // BYM2 IID component
  vector[N] u_raw;                                 // BYM2 ICAR component (raw)
  real<lower=0> sigma_b;                           // BYM2 total magnitude
  real<lower=0, upper=1> phi;                      // BYM2 mixing parameter
  vector[K] beta;                                  // covariate coefficients
}

transformed parameters {
  // Sum-to-zero constraint on ICAR component
  vector[N] u = u_raw - mean(u_raw);

  // BYM2 composite district effect (Riebler et al. 2016)
  vector[N] b = sigma_b * (
    sqrt(1 - phi) * v +
    sqrt(phi / scaling_factor) * u
  );

  // Linear predictor on the logit scale
  vector[N] theta;
  for (i in 1:N) {
    theta[i] = mu + b[i];
  }
  if (K > 0) {
    theta += X * beta;
  }
}

model {
  // Priors
  mu ~ normal(-1.85, 1);                           // moderate-burden prevalence
  sigma_b ~ exponential(4.6);                      // approximates PC prior
  phi ~ beta(3, 2);                                // mild prior toward spatial
  beta ~ normal(0, 1);                             // weakly informative on standardised covariates
  v ~ normal(0, 1);                                // BYM2 IID component

  // ICAR prior on u_raw (improper Gaussian on graph differences)
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);               // soft sum-to-zero for sampler stability

  // Fay-Herriot likelihood for observed districts only
  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  // Predicted prevalence on probability scale for all districts (observed + unobserved)
  vector<lower=0, upper=1>[N] prevalence;
  for (i in 1:N) {
    prevalence[i] = inv_logit(theta[i]);
  }

  // National-level prevalence from the intercept
  real p_national = inv_logit(mu);

  // Posterior predictive draws for observed districts (for PPC)
  vector[N_obs] y_rep;
  for (n in 1:N_obs) {
    y_rep[n] = normal_rng(theta[obs_idx[n]], psi[n]);
  }
}
