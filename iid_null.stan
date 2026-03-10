//Independent and Indentically Distributed Fay-Herroit Stan Model with no covariates

// The input data is a vector 'y' of length 'N'.
data{
  
  // dimensions of the data
  
  int<lower=0> C; //number of clusters C should be greater than 0
  int<lower=0> D; //number of districts for ghana
  
  //outcome
  array[C] int<lower=0> y; //positive cases per cluster greater or equal to 0, never negative
  array[C] int<lower=0> n; // children tested per cluster
  
  //luster to district mapping
  array[C] int<lower=1, upper=D> district; // which district does cluster c belong to
  
  // aggregation
  matrix[D, C] A;             // D x C matrix 
  vector[D] district_counts;  // number of cluster per districts
} 

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters{
  //grand mean (overall mean)
  real mu;   //overall log-odds of prevalence in ghana
  
  // random effects
  vector[D] u;         // district random effects (deviations from the grand mean IID)
  vector[C] epsilon;    // cluster deviations
  
  //variance components
  real<lower=0> sigma_u; // how much district vary around the grand mean
  real<lower=0> sigma_v; // how much cluster vary within districts
  
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model{
//Priors
mu ~ normal(0, 1);

//iid district effects - each district drawn from a normal N(0, 1)
// something like u_d <- rnorm(216, mean = 0, sd=1) in R
u ~ normal(0, 1);

// cluster overdisperions
epsilon ~ normal(0, 1);

// variance components
sigma_u ~ exponential(2);
sigma_v ~ exponential(2);

// linear predictors 
vector[C] eta = mu + sigma_u * u[district] + sigma_v * epsilon;

//likelihood
y ~ binomial_logit(n, eta);
}

generated quantities {
  vector[C] pi_c = inv_logit(
    mu + sigma_u * u[district] + sigma_v * epsilon
  );

  vector[D] pi_d = (A * pi_c) ./ district_counts;

  real pi_ghana = inv_logit(mu);

  vector[C] log_lik;
  for (c in 1:C) {
    log_lik[c] = binomial_logit_lpmf(y[c] | n[c],
                   mu + sigma_u * u[district[c]] + sigma_v * epsilon[c]);
  }
}
