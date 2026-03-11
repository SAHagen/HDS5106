// IID Model with Covariates — Model 2
// logit(pi_c) = mu + X_c'beta + sigma_u * u_d[c] + sigma_v * eps_c


data {
  int<lower=0> C;                              // number of clusters
  int<lower=0> D;                              // number of districts
  int<lower=0> P;                              // number of covariates

  array[C] int<lower=0> y;                     // malaria positives per cluster
  array[C] int<lower=0> n;                     // children tested per cluster
  array[C] int<lower=1, upper=D> district;     // district index for each cluster

  matrix[C, P] X;                              // covariate matrix: C x P
                                               // covariates must be standardized
                                               // (mean 0, sd 1) before passing

  matrix[D, C] A;                              // D x C district-cluster indicator
  vector[D] district_counts;                   // clusters per district
}

parameters {
  real mu;                        // grand mean (overall log-odds)
  vector[P] beta;                 // fixed effect coefficients (one per covariate)
  vector[D] u;                    // district random effects (IID)
  vector[C] epsilon;              // cluster overdispersion
  real<lower=0> sigma_u;          // district variance scale
  real<lower=0> sigma_v;          // cluster variance scale
}

model {
  // ---------------------------------------------------------
  // PRIORS
  // ---------------------------------------------------------

  // grand mean — N(0,1) on log-odds scale
  // logit^-1(-2) = 12%, logit^-1(2) = 88%
  // weakly informative, data will dominate
  mu ~ normal(0, 1);

  // fixed effects — N(0,1) per coefficient
  // assumes covariates are standardized (mean 0, sd 1)
  // 95% prior mass on OR in [0.14, 7.4] per 1-SD covariate change
  // global: same beta applies to all clusters regardless of district
  beta ~ normal(0, 1);

  // IID district effects
  // u_d ~ N(0,1) independently — no spatial structure
  // equivalent to: u_d <- rnorm(D, mean = 0, sd = 1) in R
  u ~ normal(0, 1);

  // cluster overdispersion
  epsilon ~ normal(0, 1);

  // variance components — Exp(2): mean = 0.5, Pr(sigma > 1) = 0.14
  // encourages shrinkage, allows large variance if data demands
  sigma_u ~ exponential(2);
  sigma_v ~ exponential(2);

  // ---------------------------------------------------------
  // LINEAR PREDICTOR (vectorized over all C clusters)
  // ---------------------------------------------------------

  // X * beta: matrix-vector product — dot product for each cluster
  // u[district]: looks up district effect for each cluster's district
  // sigma_u, sigma_v: non-centred parameterisation — separates
  //                   scale from raw effect for better sampling geometry
  vector[C] eta = mu
                  + X * beta
                  + sigma_u * u[district]
                  + sigma_v * epsilon;

  // ---------------------------------------------------------
  // LIKELIHOOD
  // ---------------------------------------------------------

  // binomial_logit applies inv_logit internally
  // numerically more stable than binomial(n, inv_logit(eta))
  y ~ binomial_logit(n, eta);
}

generated quantities {
  // cluster-level fitted prevalence
  vector[C] pi_c = inv_logit(
    mu + X * beta + sigma_u * u[district] + sigma_v * epsilon
  );

  // district-level prevalence
  // sampled districts: average over clusters
  // unsampled districts: grand mean + district random effect
  vector[D] pi_d;
  for (d in 1:D) {
    if (district_counts[d] > 0) {
      pi_d[d] = dot_product(A[d], pi_c) / district_counts[d];
    } else {
      pi_d[d] = inv_logit(mu + sigma_u * u[d]);
    }
  }

  // overall Ghana prevalence (back-transformed grand mean)
  real pi_ghana = inv_logit(mu);

  // log-likelihood for LOO model comparison
  vector[C] log_lik;
  for (c in 1:C) {
    log_lik[c] = binomial_logit_lpmf(y[c] | n[c],
                   mu
                   + dot_product(X[c], beta)
                   + sigma_u * u[district[c]]
                   + sigma_v * epsilon[c]);
  }
}
