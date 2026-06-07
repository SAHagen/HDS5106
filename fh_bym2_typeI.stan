// =============================================================================
// fh_bym2_typeI.stan
//
// Bayesian Fay-Herriot Small Area Estimation model with BYM2 spatial random
// effect, fixed year effect, and Type I (IID) space-time interaction.
//
// References:
//   - Wakefield, Okonek & Pedersen (2020): Fay-Herriot framework
//   - Riebler et al. (2016): BYM2 reparameterisation
//   - Knorr-Held (2000): Type I space-time interaction
//   - Goicoa et al. (2018): Identifiability constraints
//
// For T = 2 biomarker snapshots:
//   - The year effect uses delta[1] = 0, delta[2] estimated
//     (fixed effect, NOT random walk; appropriate for current-status biomarkers)
//   - Type I interaction: gamma_{i,t} reduces to one parameter per district
//     under the constraint gamma_{i,2} = -gamma_{i,1}
//   - Additional sum-to-zero constraint on gamma_star
// =============================================================================

data {
  int<lower=1> N;                                  // number of districts
  int<lower=1> T;                                  // number of time points (= 2)
  int<lower=1> NT;                                 // total district-time rows = N * T

  int<lower=1> N_obs;                              // number of observed rows
  array[N_obs] int<lower=1, upper=NT> obs_idx;     // indices of observed rows in panel
  array[NT] int<lower=1, upper=N> district;        // district index per row
  array[NT] int<lower=1, upper=T> time;            // time index per row

  vector[N_obs] y;                                 // logit direct estimates
  vector<lower=0>[N_obs] psi;                      // logit-scale SEs

  // Spatial adjacency graph
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  real<lower=0> scaling_factor;

  // Covariates (NT x K — can include time-varying)
  int<lower=0> K;
  matrix[NT, K] X;
}

parameters {
  real mu;                                         // national logit intercept
  vector[N] v;                                     // BYM2 IID component
  vector[N] u_raw;                                 // BYM2 ICAR component (raw)
  real<lower=0> sigma_b;                           // BYM2 total magnitude
  real<lower=0, upper=1> phi;                      // BYM2 mixing parameter
  real delta_2;                                    // year-2 fixed effect
  vector[N] gamma_raw;                             // Type I interaction (raw)
  real<lower=0> sigma_gamma;                       // interaction SD
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

  // Year effect vector (fixed effect; delta[1] = 0, delta[2] free)
  vector[T] delta;
  delta[1] = 0;
  delta[2] = delta_2;

  // Sum-to-zero constraint on space-time interaction (Goicoa et al. 2018)
  vector[N] gamma_star = gamma_raw - mean(gamma_raw);

  // Linear predictor on the logit scale
  vector[NT] theta;
  for (n in 1:NT) {
    theta[n] = mu + b[district[n]] + delta[time[n]];

    // Type I interaction with sign flip:
    //   gamma_{i,1} = +gamma_star[i],  gamma_{i,2} = -gamma_star[i]
    // Reflects the T=2 constraint gamma_{i,2} = -gamma_{i,1}
    if (time[n] == 1) {
      theta[n] += sigma_gamma * gamma_star[district[n]];
    } else {
      theta[n] -= sigma_gamma * gamma_star[district[n]];
    }
  }
  if (K > 0) {
    theta += X * beta;
  }
}

model {
  // Priors
  mu ~ normal(-1.85, 1);
  sigma_b ~ exponential(4.6);
  phi ~ beta(3, 2);
  delta_2 ~ normal(0, 5);                          // weakly informative on national trend
  sigma_gamma ~ exponential(2);                    // moderate prior on interaction SD
  beta ~ normal(0, 1);
  v ~ normal(0, 1);
  gamma_raw ~ normal(0, 1);

  // ICAR prior on u_raw
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);

  // Fay-Herriot likelihood for observed cells
  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  // Predicted prevalence on probability scale for ALL district-time cells
  vector<lower=0, upper=1>[NT] prevalence;
  for (n in 1:NT) {
    prevalence[n] = inv_logit(theta[n]);
  }

  // National prevalences for each round
  real delta_change = delta_2;
  real p_round1 = inv_logit(mu);
  real p_round2 = inv_logit(mu + delta_2);

  // Posterior predictive draws for observed cells (for PPC)
  vector[N_obs] y_rep;
  for (n in 1:N_obs) {
    y_rep[n] = normal_rng(theta[obs_idx[n]], psi[n]);
  }
}
