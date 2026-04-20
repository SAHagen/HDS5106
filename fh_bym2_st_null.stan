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

// The input data is a vector 'y' of length 'N'.
// =============================================================================
// Fay-Herriot BYM2 Space-Time Null Model
// Single country, two survey rounds
// =============================================================================
// 
// Likelihood:
//   y[i,t] | theta[i,t] ~ N(theta[i,t], psi[i,t]^2)   [observed only]
// 
// Linking model:
//   theta[i,t] = mu + b[i] + delta[t]
//
// Spatial (BYM2, Riebler et al. 2016):
//   b[i] = sigma_b * (sqrt(1-phi) * v[i] + sqrt(phi/s) * u[i])
//   v[i] ~ N(0, 1)
//   u    ~ ICAR on graph, sum(u) = 0
//
// Temporal (fixed effect, sum-to-zero for T=2):
//   delta[1] = +alpha, delta[2] = -alpha
//
// =============================================================================

data {
  int<lower=1> N;                      // number of districts
  int<lower=1> T;                      // number of survey rounds (= 2)
  int<lower=1> NT;                     // total district-year rows (= N*T)

  int<lower=1> N_obs;                  // number of observed district-years
  array[N_obs] int<lower=1, upper=NT> obs_idx;    // indices of observed rows

  array[NT] int<lower=1, upper=N> district;       // district index per row
  array[NT] int<lower=1, upper=T> time;           // time index per row

  vector[N_obs] y;                     // logit direct estimates (observed only)
  vector<lower=0>[N_obs] psi;          // known SE on logit scale (observed only)

  int<lower=0> N_edges;                // number of neighbour pairs
  array[N_edges] int<lower=1, upper=N> node1;     // edge endpoint 1
  array[N_edges] int<lower=1, upper=N> node2;     // edge endpoint 2

  real<lower=0> scaling_factor;        // BYM2 scaling factor
}

parameters {
  real mu;                             // overall intercept
  real alpha;                          // half the temporal change (delta_1)

  vector[N] v;                         // IID district effects
  vector[N] u_raw;                     // unscaled ICAR (sum-to-zero imposed below)

  real<lower=0> sigma_b;               // total district SD
  real<lower=0, upper=1> phi;          // spatial fraction of variance
  real<lower=0> sigma_delta;           // temporal SD
}

transformed parameters {
  // Sum-to-zero on u
  vector[N] u = u_raw - mean(u_raw);

  // BYM2 combined district effect
  vector[N] b = sigma_b * (sqrt(1 - phi) * v + sqrt(phi / scaling_factor) * u);

  // Temporal effects (sum-to-zero for T=2)
  vector[T] delta;
  delta[1] =  alpha;
  delta[2] = -alpha;

  // Full linear predictor for ALL district-year combinations
  vector[NT] theta;
  for (n in 1:NT) {
    theta[n] = mu + b[district[n]] + delta[time[n]];
  }
}

model {
  // =========================================================================
  // PRIORS
  // =========================================================================
  mu          ~ normal(-2, 1);                 // national baseline
  sigma_b     ~ exponential(4.6);              // PC prior on total SD
  phi         ~ beta(3, 2);                    // mild preference for spatial
  sigma_delta ~ exponential(2);                // PC prior on temporal SD
  alpha       ~ normal(0, sigma_delta);        // temporal change

  v ~ std_normal();                            // IID component

  // ICAR prior on u_raw (pairwise differences)
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  
  // Soft sum-to-zero constraint on u_raw (centred in transformed parameters)
  sum(u_raw) ~ normal(0, 0.001 * N);

  // =========================================================================
  // LIKELIHOOD (observed district-years only)
  // =========================================================================
  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  // Posterior prevalence for all district-years (observed + unobserved)
  vector<lower=0, upper=1>[NT] prevalence;
  for (n in 1:NT) {
    prevalence[n] = inv_logit(theta[n]);
  }

  // Derived summaries
  real delta_change = 2 * alpha;               // total logit change across rounds
  real p_round1 = inv_logit(mu + alpha);       // national prevalence round 1
  real p_round2 = inv_logit(mu - alpha);       // national prevalence round 2
}

