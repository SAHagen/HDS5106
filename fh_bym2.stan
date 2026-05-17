// =============================================================================
// BYM2 Fay-Herriot Model — Ghana Under-5 Malaria Prevalence (2022)
// =============================================================================
// Cross-sectional SAE with BYM2 spatial random effects
// Priors grounded in Wakefield et al. (2020), Riebler et al. (2016),
// Simpson et al. (2017), Gelman et al. (2008)
data {
  int<lower=1> N;                // Number of districts (260 for Ghana)
  int<lower=1> K;                // Number of covariates

  // Fay-Herriot inputs
  vector[N] y;                   // Direct estimates: logit(p_hat_i) (Hajek)
  vector<lower=0>[N] sigma_e;   // Known sampling SEs for each district

  // Covariate matrix (standardised)
  matrix[N, K] X;

  // Spatial adjacency structure
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;

  // BYM2 scaling factor (from INLA::inla.scale.model)
  real<lower=0> scaling_factor;

  // Prior predictive toggle: 1 = sample from priors only, 0 = full model
  int<lower=0, upper=1> prior_only;
}

parameters {
  real beta0;     // Global intercept (logit scale)
  vector[K] beta;  // Covariate coefficients
  real<lower=0> sigma;      // BYM2 total marginal SD
  real<lower=0, upper=1> phi;  // BYM2 mixing parameter (spatial fraction)
  vector[N] u;    // IID component (unscaled)
  vector[N] s;  // ICAR component (unscaled)
}

transformed parameters {
  // BYM2 convolution: combine spatial (ICAR) and non-spatial (IID)
  // Riebler et al. (2016) parameterisation:
  //   convolution_i = sigma * (sqrt(phi/scaling_factor) * s_i + sqrt(1-phi) * u_i)
  //
  // scaling_factor ensures the ICAR component has unit generalised variance,
  // so sigma is interpretable as the marginal SD of the combined effect,
  // and phi is the proportion of that variance attributable to spatial structure.

  vector[N] convolution;
  convolution = sigma * (sqrt(phi / scaling_factor) * s +
                         sqrt(1 - phi) * u);

  // Linear predictor on logit scale
  vector[N] mu;
  mu = beta0 + X * beta + convolution;
}

model {
  // Intercept — Wakefield, Okonek & Pedersen (2020)
  // N(0, sqrt(1000)) is the SUMMER/surveyPrev default.
  // Alternative: N(-2, 1) centres near Ghana's ~8.8% prevalence on logit scale.
  // Using N(-2, 1) here as it is more informative and better calibrated
  // for the prior predictive check. Defensible via Gelman et al. (2008).
  beta0 ~ normal(-2, 1);

  // Covariate effects — Gelman et al. (2008)
  // N(0, 1) on standardised covariates: a 1-SD covariate change shifts
  // log-odds by at most ~2-3, corresponding to OR in (0.05, 20).
  beta ~ normal(0, 1);

  // BYM2 total marginal SD — PC prior: P(sigma > 1) = 0.01
  // Riebler et al. (2016); Simpson et al. (2017)
  // Exponential(lambda) where lambda = -log(0.01)/1 = 4.6
  // This says: 99% prior probability that sigma < 1.
  // On the logit scale, sigma = 1 means random effects can shift
  // log-odds by ~2-3, which is already very large spatial variation.
  target += exponential_lpdf(sigma | 4.6);

  // BYM2 mixing parameter — Riebler et al. (2016)
  // PC prior approximated as Beta(0.5, 0.5).
  // P(phi > 0.5) = 2/3: gently favours spatial structure,
  // which is reasonable for disease mapping where neighbouring
  // districts share risk factors.
  phi ~ beta(0.5, 0.5);

  // ICAR prior on spatial component s
  // The pairwise difference penalty induces spatial smoothing:
  // neighbouring districts get similar s values.
  // This is the intrinsic CAR (Besag, 1974; Besag et al., 1991).
  target += -0.5 * dot_self(s[node1] - s[node2]);

  // Soft sum-to-zero constraint on s for identifiability
  // (the ICAR is rank-deficient; this anchors the level)
  sum(s) ~ normal(0, 0.001 * N);

  // IID component: standard normal (unscaled; sigma handles the scale)
  u ~ normal(0, 1);


  // ===================================================================
  // LIKELIHOOD — only included when prior_only == 0
  // ===================================================================
  // Fay-Herriot: the direct estimate y_i is normally distributed
  // around the true logit-prevalence mu_i, with known SE sigma_e_i.
  //
  // y_i ~ Normal(mu_i, sigma_e_i)
  //
  // This is the area-level SAE model (Fay & Herriot, 1979),
  // recommended by Wakefield et al. (2020) as the most reliable
  // model choice for DHS prevalence mapping because it acknowledges
  // the survey design through the weighted estimate and its variance.

  if (prior_only == 0) {
    y ~ normal(mu, sigma_e);
  }
}

generated quantities {
  // District-level prevalence on probability scale
  vector[N] prev;
  prev = inv_logit(mu);

  // Simulated direct estimates (prior or posterior predictive)
  // These are what the data WOULD look like under the current
  // parameter values. Compare y_rep to y for model checking.
  vector[N] y_rep;
  for (i in 1:N) {
    y_rep[i] = normal_rng(mu[i], sigma_e[i]);
  }

  // Log-likelihood (for LOO-CV, only meaningful when prior_only == 0)
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | mu[i], sigma_e[i]);
  }
}
