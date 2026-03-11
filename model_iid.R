library(sf)
library(terra)
library(geodata)
library(exactextractr)
library(spdep)
library(tmap)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(haven)
library(survey)
library(srvyr)
library(malariaAtlas)
library(rstan)
library(loo)
library(Matrix)

# DHS cluster GPS points
dhs_clusters <- read_sf("dhs_gh/GHGE8AFL.shp")

# Ghana district polygons from GADM (admin level 2)
ghana_districts <- gadm(country = "GHA", level = 2, path = "dhs_gh/") %>%
  st_as_sf() %>%
  st_transform(crs = 4326) %>%
  rename(district_name = NAME_2, district_code = GID_2) 

# Standardize the cluster GPS data and make sure CRS matches districts
gps <- dhs_clusters %>%
  rename(cluster_id = DHSCLUST) %>%
  st_transform(crs = 4326)

# both layers should be in the same CRS before joining
stopifnot(st_crs(gps) == st_crs(ghana_districts)) 

#Assigning each DHS cluster to its nearest district using  the nearest feature
gps_joined <- st_join(gps, ghana_districts, join = st_nearest_feature) %>%
  select(cluster_id, district_name, district_code) %>%
  st_drop_geometry()


# Building a district lookup table (name + code → numeric ID for Stan)
district_lookup <- ghana_districts %>%
  st_drop_geometry() %>%
  select(district_name, district_code) %>%
  distinct() %>%
  arrange(district_name) %>%
  mutate(district_id = row_number())

# Attaching numeric district IDs to the cluster-district mapping
gps_joined <- gps_joined %>%
  left_join(district_lookup, by = c("district_name", "district_code"))

#Loading and cleaning DHS household recode file
# This file has one row per person; we keep children under 5 (60 months) with RDT results
dhs_hm <- read_dta("dhs_gh/GHPR8CFL.dta")

pr <- dhs_hm %>%
  select(
    cluster_id = hv001,
    hh_id = hv002,
    line_id = hvidx,
    weight_raw  = hv005,
    stratum  = hv023,
    age_months= hc1,
    malaria_rdt = hml32,
    net_use  = hml12,
    wealth = hv270,
    urban_rural  = hv025
  ) %>%
  filter(age_months < 60) %>%           # keep children under 5 only
  filter(malaria_rdt %in% c(0, 1)) %>%  # keep only valid RDT results (drops missing/refused)
  mutate(
    weight= weight_raw / 1e6,       # DHS weights are scaled by 1,000,000
    malaria = as.integer(malaria_rdt), 
    cluster_id = as.integer(cluster_id),
    stratum = as.integer(stratum)
  ) %>%
  select(-weight_raw, -malaria_rdt)


# Mergeing individual data with district assignments
my_df <- pr %>%
  left_join(gps_joined, by = "cluster_id")

# Survey design (DHS uses stratified two-stage cluster sampling)
# Stage 1: clusters within strata (region × urban/rural)
# Stage 2: households within clusters
# Weights are inverse-probability weights provided by DHS
dhs_design <- svydesign(
  ids  = ~cluster_id,
  strata  = ~stratum,
  weights = ~weight,
  data= my_df,
  nest = TRUE
)

# National prevalence estimate, ~8.6% to match official Ghana 2022 DHS report
overall_prev <- svymean(~malaria, design = dhs_design)
print(overall_prev)


# Aggregating to cluster level for Stan
cluster_df <- my_df %>%
  group_by(cluster_id, district_id, district_name, stratum) %>%
  summarise(
    n_raw = n(),
    y_raw = sum(malaria),
    w_sum = sum(weight),
    yw_sum = sum(weight * malaria),
    .groups = "drop"
  ) %>%
  mutate(
    prev_w = yw_sum / w_sum,  # weighted prevalence per cluster
    n = as.integer(n_raw),
    y = as.integer(y_raw),
    y = pmax(0L, pmin(y, n))# y can't exceed n or go below 0
  ) %>%
  filter(n > 0)

#Building Stan data list
C <- nrow(cluster_df) # total number of clusters
D <- nrow(district_lookup)# total number of districts

# A is a D × C indicator matrix: A[d, c] = 1 if cluster c belongs to district d
# Used in Stan to aggregate cluster-level estimates up to district level
A <- sparseMatrix(
  i = cluster_df$district_id,
  j = 1:C,
  x = 1,
  dims = c(D, C)
) %>% as.matrix()

district_counts <- rowSums(A)  # how many clusters each district has

stan_data <- list(
  C = C,
  D = D,
  y = cluster_df$y,
  n  = cluster_df$n,
  district  = cluster_df$district_id,
  A  = A,
  district_counts = district_counts
)

#Fitting
options(mc.cores = parallel::detectCores())

fit_iid <- stan(
  file  = "iid_null.stan",
  data = stan_data,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

# summary of key parameters from the fitted model
print(fit_iid,
      pars  = c("mu", "sigma_u", "sigma_v", "pi_ghana"),
      probs = c(0.025, 0.5, 0.975))

#Extracting posterior estimates 
posterior <- extract(fit_iid)
pi_d_samples <- posterior$pi_d # pi_d is the posterior distribution of malaria prevalence for each district

pi_d_summary <- tibble(
  district_id = 1:D,
  mean  = colMeans(pi_d_samples),
  median = apply(pi_d_samples, 2, median),
  lower = apply(pi_d_samples, 2, quantile, 0.025),
  upper = apply(pi_d_samples, 2, quantile, 0.975),
  ci_width = upper - lower   # uncertainty spread
)

# Attaching estimates to Ghana map
ghana_map_estimates <- ghana_districts %>%
  left_join(district_lookup, by = c("district_name", "district_code")) %>%
  left_join(pi_d_summary,    by = "district_id")

#Raw prevalence by district for comparison
district_raw <- cluster_df %>%
  group_by(district_id, district_name) %>%
  summarise(
    raw_prev  = sum(y) / sum(n),
    n_children = sum(n),
    n_clusters = n(),
    .groups = "drop"
  )

ghana_map_raw <- ghana_districts %>%
  left_join(district_lookup, by = c("district_name", "district_code")) %>%
  left_join(district_raw,    by = c("district_id", "district_name"))

#Maps: raw vs model estimates
map_raw <- tm_shape(ghana_map_raw) +
  tm_polygons(
    col      = "raw_prev",
    palette  = "-magma",
    title    = "Raw",
    alpha    = 0.8,
    breaks   = seq(0, 0.4, by = 0.05),
    colorNA  = "grey85",
    textNA   = "No data",
    popup.vars = c("District" = "district_name", "Raw" = "raw_prev")
  )

map_model <- tm_shape(ghana_map_estimates) +
  tm_polygons(
    col      = "mean",
    palette  = "-magma",
    title    = "Model",
    alpha    = 0.8,
    breaks   = seq(0, 0.4, by = 0.05),
    colorNA  = "grey85",
    textNA   = "No data",
    popup.vars = c("District" = "district_name", "Posterior" = "mean")
  )

tmap_arrange(map_raw, map_model, ncol = 2)

############ 


# Model 2 with covariates 

covariates_data <- read_csv("dhs_gh/GHGC8AFL.csv") %>% 
  rename(cluster_id = DHSCLUST) %>% 
  select(
    cluster_id,
    rainfall   = Rainfall_2020,
    elevation= Elevation,
    nightlights = Nightlights_Composite,
    itn_coverage= ITN_Coverage_2020,
    temp_range = Diurnal_Temperature_Range_2020,
    travel_time = Travel_Times)

# joining the covatiates to cluster _df created before

cluster_df_cov <- cluster_df %>% 
  left_join(covariates_data, by ="cluster_id")

# standardizing the covariates
cluster_df_cov <- cluster_df_cov %>%
  mutate(across(c(rainfall, elevation, nightlights,
                  itn_coverage, temp_range, travel_time),
                ~as.numeric(scale(.)),
                .names = "{.col}_std"))

# creating the covaritates matrix

X <- cluster_df_cov %>% 
  select(rainfall_std, elevation_std, nightlights_std,
         itn_coverage_std, temp_range_std, travel_time) %>%
  as.matrix()
  

# appending the covariate matrix to the stan model
stan_data_cov <- stan_data
stan_data_cov$X <- X # teh covatiates matrix defined above
stan_data_cov$P <- ncol(X) # number of covatiates

# fitting the second model with covariates

fit_iid_cov <- stan(
  file  = "iid_covariates.stan",
  data = stan_data_cov,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

print(fit_iid_cov,
      pars  = c("mu", "sigma_u", "sigma_v", "pi_ghana"),
      probs = c(0.025, 0.5, 0.975))

print(fit_iid_cov,
      pars  = "beta",
      probs = c(0.025, 0.5, 0.975))
