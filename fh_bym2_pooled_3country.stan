// Pooled three-country BYM2 Fay-Herriot space-time model
// with country-specific random slopes and cross-border ICAR
//
// Features:
//   - Single BYM2 over 724-district graph (cross-border smoothing)
//   - Country fixed effects (intercept shift)
//   - Country-specific time effects delta_c (CIV goes up, others down)
//   - Hierarchical random slopes: beta_ck ~ N(beta_bar_k, tau_k)
//   - Knorr-Held Type I IID interaction
//
// Author: Godwill Zulu | May 2026

data {
  int<lower=1> N;                    // total districts (724)
  int<lower=1> T;                    // time points (2)
  int<lower=1> NT;                   // N * T (1448)
  int<lower=1> C;                    // number of countries (3)
  int<lower=1> K;                    // number of covariates (7)
  
  // Observation structure
  int<lower=1> N_obs;                // observed district-time pairs
  array[N_obs] int<lower=1, upper=NT> obs_idx;  // which rows are observed
  vector[N_obs] y;                   // logit direct estimates (observed only)
  vector<lower=0>[N_obs] psi;        // known sampling SE (observed only)
  
  // Panel structure
  array[NT] int<lower=1, upper=N> district;  // district index for each row
  array[NT] int<lower=1, upper=T> time;      // time index for each row
  array[N] int<lower=1, upper=C> country;    // country index for each district
  
  // Covariates
  matrix[NT, K] X;                   // covariate matrix (all district-times)
  
  // Spatial structure (combined 724-district graph)
  int<lower=1> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  real<lower=0> scaling_factor;
}

parameters {
  // Global intercept
  real mu;
  
  // Country effects (C-1 free, corner constraint: alpha[1] = 0)
  vector[C - 1] alpha_raw;
  
  // Country-specific time effects (delta_{c,1} = 0, estimate delta_{c,2})
  vector[C] delta_2;
  
  // BYM2 spatial (over full 724-district graph)
  vector[N] v;                       // IID component
  vector[N] u_raw;                   // ICAR component (raw)
  real<lower=0> sigma_b;             // total spatial SD
  real<lower=0, upper=1> phi;        // mixing parameter
  
  // Hierarchical random slopes
  matrix[C, K] beta_raw;            // country-specific slopes (raw)
  vector[K] beta_bar;               // global mean slopes
  vector<lower=0>[K] tau;           // between-country SD per covariate
  
  // Knorr-Held Type I interaction
  vector[NT] gamma_raw;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  // Country intercepts (corner constraint)
  vector[C] alpha;
  alpha[1] = 0;
  alpha[2:C] = alpha_raw;
  
  // Country-specific slopes
  matrix[C, K] beta;
  for (k in 1:K)
    for (c in 1:C)
      beta[c, k] = beta_bar[k] + tau[k] * beta_raw[c, k];
  
  // BYM2 spatial effects
  vector[N] u = u_raw - mean(u_raw);
  vector[N] b = sigma_b * (sqrt(1 - phi) * v + sqrt(phi / scaling_factor) * u);
  
  // Interaction effects (soft sum-to-zero)
  vector[NT] gamma_star = sigma_gamma * gamma_raw;
  
  // Linear predictor for ALL district-time pairs
  vector[NT] theta;
  for (n in 1:NT) {
    int d = district[n];
    int t = time[n];
    int c = country[d];
    
    // Country-specific covariate effects
    real xb = 0;
    for (k in 1:K)
      xb += X[n, k] * beta[c, k];
    
    theta[n] = mu + alpha[c] + b[d]
               + (t == 2 ? delta_2[c] : 0)
               + gamma_star[n]
               + xb;
  }
}

model {
  // ── Priors ──────────────────────────────────────────────────────
  mu ~ normal(-1.85, 1.5);           // slightly wider for multi-country
  alpha_raw ~ normal(0, 1.5);        // country intercept shifts
  delta_2 ~ normal(0, 2);            // country-specific time effects
  
  // BYM2
  sigma_b ~ exponential(4.6);        // PC prior
  phi ~ beta(3, 2);
  v ~ std_normal();
  target += -0.5 * dot_self(u_raw[node1] - u_raw[node2]);
  sum(u_raw) ~ normal(0, 0.001 * N);
  
  // Hierarchical slopes
  beta_bar ~ normal(0, 1);           // global mean effects
  tau ~ exponential(2);              // shrinkage toward shared effects
  to_vector(beta_raw) ~ std_normal();  // non-centred parameterisation
  
  // Interaction
  sigma_gamma ~ exponential(2);
  gamma_raw ~ std_normal();
  
  // Soft sum-to-zero for gamma (per district and per time)
  // Not strictly necessary with Type I but helps identification
  
  // ── Likelihood (observed only) ─────────────────────────────────
  y ~ normal(theta[obs_idx], psi);
}

generated quantities {
  // Log-likelihood for LOO-CV
  vector[N_obs] log_lik;
  for (i in 1:N_obs)
    log_lik[i] = normal_lpdf(y[i] | theta[obs_idx[i]], psi[i]);
  
  // Prevalence on probability scale (all districts, all times)
  vector[NT] prevalence;
  for (n in 1:NT)
    prevalence[n] = inv_logit(theta[n]);
  
  // National prevalence per country per round
  // (average prevalence across districts within each country)
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
  
  // Country-specific temporal change (odds ratio)
  vector[C] odds_ratio;
  for (c in 1:C)
    odds_ratio[c] = exp(delta_2[c]);
}
