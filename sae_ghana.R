library(tidyverse)    
library(sf)           # Simple features for spatial data
library(terra)        # Modern raster handling
library(exactextractr)# Extract raster values by polygon
library(spdep)        # Spatial dependence and neighbors
library(haven)        # Read Stata (.dta) files
library(survey)       # Survey analysis with complex designs
library(srvyr)        # Tidyverse-friendly survey analysis
library(malariaAtlas) # Download MAP rasters
library(rstan)        # Interface to Stan
library(loo)          # Model comparison
library(tmap)         # Thematic maps
library(ggplot2)      # Plots
library(patchwork)    # Combine plots

# Stan configuration
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Data Loading and Prepration from DHS and WHo boundaries 
who_boundaries <- read_sf("glob_polygons/GLOBAL_ADM2.shp")
dhs_clusters <- read_sf("dhs_gh/GHGE8AFL.shp") # HS Clusters
dhs_children <- read_dta("dhs_gh/GHKR8CFL.dta") # KR file for ACT Uptake
dhs_hm <- read_dta("dhs_gh/GHPR8CFL.dta") # PR file for ITN and Malaria Results

ghana_districts <- who_boundaries %>%
  filter(ADM0_NAME == "GHANA") %>% # filtering for ghana and transfomrin
  st_transform(crs = 4326)

gps <- dhs_clusters %>% 
  rename(cluster_id = DHSCLUST) # renaming cluster id for easy joins and use

gps_joined <- st_join(gps, ghana_districts, join = st_nearest_feature) # joining bounderies and dhs clusters 

# KR file (Children's data) for ACT intervention
kr <- dhs_children %>% 
  select(
    cluster_id = v001,
    hh_id = v002,
    line_id = b16,
    fever = h22,
    took_act = h32z
  ) %>% 
  mutate(
    act_binary = case_when(
      took_act == 1 ~ 1,
      took_act == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )
# PR file (Household member data)ITN, and Malaria Outcome
pr <- dhs_hm %>% 
  select(
    cluster_id = hv001,
    hh_id = hv002,
    line_id = hvidx,
    weight_raw = hv005,
    stratum = hv023,
    age_months = hc1,
    malaria_results = hml35,
    net_use = hml12,
    wealth = hv270,
    urban_rural = hv025
  ) %>% 
  filter(age_months < 60) %>%  # Children under 5
  mutate(
    weight = weight_raw / 1000000,
    outcome_binary = case_when(
      malaria_results == 1 ~ 1,
      malaria_results == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    itn_binary = case_when(
      net_use == 1 ~ 1,
      net_use == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

#MJoining datasets
my_data_dhs <- kr %>%
  left_join(pr, by = c("cluster_id", "hh_id", "line_id"))

my_df <- my_data_dhs %>%
  left_join(st_drop_geometry(gps_joined), by = "cluster_id")

# Calculating direct estimates 
# survey design
survey_design <- my_df %>%
  filter(!is.na(weight)) %>%
  as_survey_design(
    ids = cluster_id,
    strata = stratum,
    weights = weight
  )

#district level estimates
district_estimates <- survey_design %>%
  group_by(ADM2_NAME) %>% 
  summarise(
    malaria_prevalence = survey_mean(outcome_binary, vartype = c("ci", "se"), na.rm = TRUE),
    itn_coverage = survey_mean(itn_binary, vartype = c("ci", "se"), na.rm = TRUE),
    n_obs = n(),
    n_positive = sum(outcome_binary, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    # Coefficient of Variation
    malaria_cv = malaria_prevalence_se / malaria_prevalence,
    # Reliability classification
    malaria_est_rel = case_when(
      malaria_cv < 0.20 ~ "Reliable",
      malaria_cv >= 0.20 & malaria_cv < 0.30 ~ "Caution",
      malaria_cv >= 0.30 ~ "Unreliable",
      is.na(malaria_cv) | is.infinite(malaria_cv) ~ "Undefined"
    )
  )

#Gettign covariates data from MAP
ghana_shp <- getShp(ISO = "GHA", admin_level = "admin0")

#PfPR raster
pfpr_raster <- tryCatch({
  getRaster(
    dataset_id = "Malaria__202206_Global_Pf_Parasite_Rate",
    shp = ghana_shp,
    year = 2020
  )
})

#Accessibility raster
access_raster <- tryCatch({
  getRaster(
    dataset_id = "Accessibility__202001_Global_Motorized_Travel_Time_to_Healthcare",
    shp = ghana_shp
  )
})

#Extract district-level covariate summaries 
districts_for_extract <- ghana_districts %>%
  filter(ADM2_NAME %in% district_estimates$ADM2_NAME)

# Extract PfPR
if (!is.null(pfpr_raster)) {
  pfpr_values <- exact_extract(pfpr_raster, districts_for_extract, fun = "mean")
} else {
  pfpr_values <- rep(NA, nrow(districts_for_extract))
}

# Extract Accessibility
if (!is.null(access_raster)) {
  access_values <- exact_extract(access_raster, districts_for_extract, fun = "mean")
} else {
  access_values <- rep(NA, nrow(districts_for_extract))
}

# Creating covariates dataframe
district_covariates <- data.frame(
  ADM2_NAME = districts_for_extract$ADM2_NAME,
  pfpr_map = pfpr_values,
  accessibility = access_values
)

#Merge covariates with direct estimates
sae_data <- district_estimates %>%
  left_join(district_covariates, by = "ADM2_NAME") %>%
  filter(!is.na(malaria_prevalence) & !is.na(malaria_prevalence_se)) %>%
  filter(malaria_prevalence_se > 0)

#Data Preparation for stan

#cleaning and tranforming 
sae_data_clean <- sae_data %>%
  filter(
    !is.na(malaria_prevalence),
    !is.na(malaria_prevalence_se),
    malaria_prevalence_se > 0
  ) %>%
  mutate(
    # Bounding prevalence away from 0 and 1
    p_bounded = pmax(0.02, pmin(0.98, malaria_prevalence)),
    
    # Logit transformation
    logit_prev = log(p_bounded / (1 - p_bounded)),
    
    # Delta method variance transformation with capping
    logit_var = (malaria_prevalence_se^2) / 
                (pmax(0.02, malaria_prevalence)^2 * 
                 pmax(0.02, 1 - malaria_prevalence)^2),
    
    # Caping SE between 0.1 and 2.0 to avoid extremes
    logit_se = pmax(0.1, pmin(sqrt(logit_var), 2.0)),
    
    # Standardizing covariates
    itn_std = as.numeric(scale(itn_coverage)),
    pfpr_std = as.numeric(scale(pfpr_map)),
    access_std = as.numeric(scale(accessibility))
  ) %>%
  # Handling any remaining NAs in standardized covariates (they making it chard for the chains to merge)
  mutate(
    itn_std = if_else(is.na(itn_std), 0, itn_std),
    pfpr_std = if_else(is.na(pfpr_std), 0, pfpr_std),
    access_std = if_else(is.na(access_std), 0, access_std)
  ) %>%
  filter(is.finite(logit_prev), is.finite(logit_se))

# Creating spatial adjacency matrix 
map_districts <- ghana_districts %>%
  filter(ADM2_NAME %in% sae_data_clean$ADM2_NAME) %>%
  arrange(match(ADM2_NAME, sae_data_clean$ADM2_NAME))

#need to verify alignment if need be!!!!!
# Creating neighborhood structure for the Spatial ICAR model
nb <- poly2nb(map_districts, queen = TRUE)
W <- nb2mat(nb, style = "B", zero.policy = TRUE)

# Creating edge list for Stan
edges <- which(W == 1, arr.ind = TRUE)
edges <- edges[edges[,1] < edges[,2], ]  # Upper triangle only

n_districts <- nrow(sae_data_clean)
n_edges <- nrow(edges)
num_neighbors <- rowSums(W)


#Stan model specifications
#staring off with a simple abse modelw ith no covariates to test for convergence and data appropriateness
stan_model_iid_simple <- "
data {
  int<lower=1> N;              // Number of districts
  vector[N] y;                 // Direct estimates (logit scale)
  vector<lower=0>[N] psi;      // Known sampling SEs
}

parameters {
  real mu;                     // Overall mean (intercept)
  real<lower=0> sigma_u;       // SD of random effects
  vector[N] u_raw;             // Raw random effects (non-centered)
}

transformed parameters {
  vector[N] theta;
  // Non-centered parameterization: u = sigma_u * u_raw
  theta = mu + sigma_u * u_raw;
}

model {
  // Priors
  mu ~ normal(-1.5, 0.5);      // ~18% baseline prevalence
  sigma_u ~ exponential(2);    // Encourages smaller values
  u_raw ~ std_normal();        // Non-centered
  
  // Likelihood
  y ~ normal(theta, psi);
}

generated quantities {
  vector[N] prevalence;
  vector[N] log_lik;
  
  for (i in 1:N) {
    prevalence[i] = inv_logit(theta[i]);
    log_lik[i] = normal_lpdf(y[i] | theta[i], psi[i]);
  }
}
"

# IID with Covariates
stan_model_iid_cov <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] psi;
  int<lower=1> P;              // Number of covariates (including intercept)
  matrix[N, P] X;              // Design matrix
}

parameters {
  vector[P] beta;              // Regression coefficients
  real<lower=0> sigma_u;
  vector[N] u_raw;
}

transformed parameters {
  vector[N] theta;
  theta = X * beta + sigma_u * u_raw;
}

model {
  beta ~ normal(0, 1);
  sigma_u ~ exponential(2);
  u_raw ~ std_normal();
  y ~ normal(theta, psi);
}

generated quantities {
  vector[N] prevalence;
  vector[N] log_lik;
  
  for (i in 1:N) {
    prevalence[i] = inv_logit(theta[i]);
    log_lik[i] = normal_lpdf(y[i] | theta[i], psi[i]);
  }
}
"
#Spatial ICAR
stan_model_spatial <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] psi;
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
}

parameters {
  real mu;
  real<lower=0> sigma_phi;
  vector[N] phi_raw;
}

transformed parameters {
  vector[N] phi;
  vector[N] theta;
  
  // Sum-to-zero constraint
  phi = phi_raw - mean(phi_raw);
  theta = mu + sigma_phi * phi;
}

model {
  // Priors
  mu ~ normal(-1.5, 0.5);
  sigma_phi ~ exponential(2);
  
  // ICAR prior: penalize differences between neighbors
  target += -0.5 * dot_self(phi_raw[node1] - phi_raw[node2]);
  
  // Soft sum-to-zero
  sum(phi_raw) ~ normal(0, 0.001 * N);
  
  // Likelihood
  y ~ normal(theta, psi);
}

generated quantities {
  vector[N] prevalence;
  vector[N] log_lik;
  
  for (i in 1:N) {
    prevalence[i] = inv_logit(theta[i]);
    log_lik[i] = normal_lpdf(y[i] | theta[i], psi[i]);
  }
}
"

# Data for simple IID model
stan_data_simple <- list(
  N = n_districts,
  y = sae_data_clean$logit_prev,
  psi = sae_data_clean$logit_se
)

# Data for IID with covariates
X_matrix <- model.matrix(~ 1 + itn_std + pfpr_std + access_std, 
                          data = sae_data_clean)

stan_data_cov <- list(
  N = n_districts,
  y = sae_data_clean$logit_prev,
  psi = sae_data_clean$logit_se,
  P = ncol(X_matrix),
  X = X_matrix
)

# Data for spatial model
stan_data_spatial <- list(
  N = n_districts,
  y = sae_data_clean$logit_prev,
  psi = sae_data_clean$logit_se,
  N_edges = n_edges,
  node1 = edges[, 1],
  node2 = edges[, 2]
)

#Compiling models
model_iid_simple <- stan_model(model_code = stan_model_iid_simple)
model_iid_cov <- stan_model(model_code = stan_model_iid_cov)
model_spatial <- stan_model(model_code = stan_model_spatial)

#fitting model 1
fit_iid_simple <- sampling(
  model_iid_simple,
  data = stan_data_simple,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 12345,
  control = list(adapt_delta = 0.95)
)

print(fit_iid_simple, pars = c("mu", "sigma_u"))

# fitting model 2
fit_iid_cov <- sampling(
  model_iid_cov,
  data = stan_data_cov,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 12345,
  control = list(adapt_delta = 0.95)
)

print(fit_iid_cov, pars = c("beta", "sigma_u"))

# model 3
fit_spatial <- sampling(
  model_spatial,
  data = stan_data_spatial,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 12345,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit_spatial, pars = c("mu", "sigma_phi"))

# Model comparison using the Leave One Out Cross Validation ( leave one district out and try to predeicts its prevalence using the 136 remaining districts)
loo_simple <- loo(fit_iid_simple)
loo_cov <- loo(fit_iid_cov)
loo_spatial <- loo(fit_spatial)

print(loo_compare(loo_simple, loo_cov, loo_spatial))

# results comparison 
extract_prevalence <- function(fit, model_name) {
  prev_samples <- extract(fit, pars = "prevalence")$prevalence
  
  data.frame(
    district_id = 1:ncol(prev_samples),
    sae_prev = colMeans(prev_samples),
    sae_se = apply(prev_samples, 2, sd),
    sae_low = apply(prev_samples, 2, quantile, 0.025),
    sae_upp = apply(prev_samples, 2, quantile, 0.975),
    model = model_name
  )
}

# ---- 9.2 Extract from all models ----
results_simple <- extract_prevalence(fit_iid_simple, "IID Simple")
results_cov <- extract_prevalence(fit_iid_cov, "IID + Covariates")
results_spatial <- extract_prevalence(fit_spatial, "Spatial ICAR")

# comparison dataframe
comparison <- sae_data_clean %>%
  mutate(district_id = row_number()) %>%
  left_join(results_simple %>% select(district_id, 
                                       iid_prev = sae_prev, 
                                       iid_se = sae_se), 
            by = "district_id") %>%
  left_join(results_spatial %>% select(district_id, 
                                        spatial_prev = sae_prev, 
                                        spatial_se = sae_se), 
            by = "district_id") %>%
  mutate(
    # CVs
    direct_cv = malaria_prevalence_se / malaria_prevalence,
    iid_cv = iid_se / iid_prev,
    spatial_cv = spatial_se / spatial_prev,
    
    # Reliability
    direct_reliable = direct_cv < 0.20,
    iid_reliable = iid_cv < 0.20,
    spatial_reliable = spatial_cv < 0.20
  )

#Plots
# SE Comparison Scatter Plots
p1 <- ggplot(comparison, aes(x = malaria_prevalence_se, y = iid_se)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(size = n_obs), alpha = 0.6, color = "steelblue") +
  labs(
    title = "A. Direct vs IID SAE",
    x = "Direct SE",
    y = "IID SAE SE",
    size = "Sample Size"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

p2 <- ggplot(comparison, aes(x = malaria_prevalence_se, y = spatial_se)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(size = n_obs), alpha = 0.6, color = "darkgreen") +
  labs(
    title = "B. Direct vs Spatial SAE",
    x = "Direct SE",
    y = "Spatial SAE SE",
    size = "Sample Size"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

se_comparison_plot <- p1 + p2
ggsave("SE_comparison.png", se_comparison_plot, width = 12, height = 5, dpi = 300)


#Sample Size vs SE Plot
p3 <- ggplot(comparison, aes(x = n_obs, y = malaria_prevalence_se)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(NA, 0.8)) +
  labs(
    title = "High Uncertainty in Small Samples",
    x = "Sample Size (Number of Children)",
    y = "Standard Error (Uncertainty)"
  ) +
  theme_minimal()

p4 <- ggplot(district_estimates, aes(x = n_obs, y = malaria_prevalence)) +
  geom_point(alpha = 0.6, color = "#2c7fb8") +
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 0.5) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(district_estimates$n_obs, na.rm=TRUE), y = 0.105, 
           label = "Unreliable Threshold", hjust = 1, size = 3) +
  labs(
    title = "High Uncertainty in Small Samples",
    x = "Sample Size (Number of Children)",
    y = "Standard Error (Uncertainty)"
  ) +
  theme_minimal()

se_final_comparison_plot <- p4 + p3
ggsave("se_fin_comparison.png", se_final_comparison_plot, width = 8, height = 6, dpi = 300)

# Prepare map data
map_data <- map_districts %>%
  left_join(
    comparison %>% select(ADM2_NAME, 
                          direct_prev = malaria_prevalence,
                          iid_prev, spatial_prev),
    by = "ADM2_NAME"
  )

# Common scale
prev_max <- max(c(map_data$direct_prev, map_data$iid_prev, map_data$spatial_prev), 
                na.rm = TRUE)

#maps
map_direct <- tm_shape(map_data) +
  tm_polygons(
    fill = "direct_prev",
    fill.scale = tm_scale_continuous(values = "brewer.yl_or_rd", 
                                      limits = c(0, prev_max)),
    fill.legend = tm_legend(title = "Prevalence")
  ) +
  tm_title("A. Direct Estimates") +
  tm_layout(frame = FALSE)

map_iid <- tm_shape(map_data) +
  tm_polygons(
    fill = "iid_prev",
    fill.scale = tm_scale_continuous(values = "brewer.yl_or_rd", 
                                      limits = c(0, prev_max)),
    fill.legend = tm_legend(title = "Prevalence")
  ) +
  tm_title("B. IID SAE Estimates") +
  tm_layout(frame = FALSE)

map_spatial <- tm_shape(map_data) +
  tm_polygons(
    fill = "spatial_prev",
    fill.scale = tm_scale_continuous(values = "brewer.yl_or_rd", 
                                      limits = c(0, prev_max)),
    fill.legend = tm_legend(title = "Prevalence")
  ) +
  tm_title("C. Spatial SAE Estimates") +
  tm_layout(frame = FALSE)

combined_maps <- tmap_arrange(map_direct, map_iid, map_spatial, ncol = 3)
tmap_save(combined_maps, "prevalence_comparison_maps.png", 
          width = 15, height = 5, dpi = 300)



