# Multilevel Hierarchical Models for Malaria Small Area Estmation in Stan Model definations

# Model 1, IID with no ovariates 
library(tidyverse)    
library(sf)           
library(terra)        
library(exactextractr)
library(spdep)        
library(haven)        
library(survey)       
library(srvyr)        
library(malariaAtlas) 
library(rstan)      
library(loo)          
library(tmap)       
library(ggplot2)      
library(patchwork)    
library(geodata)


dhs_clusters <- read_sf("dhs_gh/GHGE8AFL.shp")

# taking ghana districts polygones from GADM
ghana_districts <- gadm(
  country = "GHA",
  level   = 2,
  path    = "dhs_gh/"
) %>%
  st_as_sf() %>%
  st_transform(crs = 4326)

gps <- dhs_clusters %>% 
  rename(cluster_id = DHSCLUST) %>%  # renaming cluster id for easy joins and use
  st_transform(crs = 4326)

# checking if the coordinates align 
st_crs(gps) == st_crs(ghana_districts)


gps_joined <- st_join(
  gps, ghana_districts, join = st_nearest_feature)  %>%
  select(
    cluster_id,
    district_name = NAME_2,    # district name from GADM
    district_code = GID_2      # unique GADM district identifier
  ) %>%
  st_drop_geometry() 


# district look up
district_lookup <- gps_joined %>% 
  distinct(district_name, district_code) %>% 
  arrange(district_name) %>% 
  mutate(district_id = row_number())
  
gps_joined <-  gps_joined %>% 
  left_join(district_lookup, by = c("district_name", "district_code"))



# laoding in and cleaning the household file from dhs with malaria test results 
dhs_hm <- read_dta("dhs_gh/GHPR8CFL.dta")
nrow(dhs_hm) # 69684 total observations


pr <- dhs_hm %>% 
  select(
    cluster_id = hv001,
    hh_id = hv002,
    line_id = hvidx,
    weight_raw = hv005,
    stratum = hv023,
    age_months = hc1,
    malaria_rdt = hml32,
    net_use = hml12,
    wealth = hv270,
    urban_rural = hv025
  ) %>% 
  filter(age_months < 60) %>%  # Children under 5
  filter(malaria_rdt %in% c(0, 1)) %>%
  mutate(
    weight = weight_raw / 1000000,
    malaria = case_when(
      malaria_rdt == 1 ~ 1L,
      malaria_rdt == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    cluster_id = as.integer(cluster_id),
    stratum = as.integer(stratum)
  ) %>%
  select(-weight_raw, -malaria_rdt) 


my_df <- pr %>%
  left_join(st_drop_geometry(gps_joined), by = "cluster_id")

# survey design and cluster level aggregation
# dhs uses stratified two stage cluster sampling, stage 1 is the clusters sampled within strata
# stage 2 is households sampled within clusters
# they defined ids -> PSU, strata -> sampling strata ( region, rural or urban)
# weights -> inverse probabilit weights, 
dhs_design <- svydesign(
  ids = ~cluster_id,
  strata= ~stratum,
  weights= ~weight,
  data = my_df,
  nest= TRUE
)

# overall malaria ghana overall prevealnce based on these dhs raw estimates 
overall_prev <- svymean(~ malaria, design = dhs_design)

overall_prev ## 8.6% as reported by ghana statistics in 2023

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
    prev_w = yw_sum / w_sum,
    n = as.integer(n_raw),
    y = as.integer(y_raw),
    y = pmax(0L, pmin(y, n))
  ) %>%
  filter(n > 0)


## Stan data set up
# dimensions
C <- nrow(cluster_df)
D <- n_distinct(cluster_df$district_id)

# district indicator matrix D x C
A_sparse <- sparseMatrix(
  i = cluster_df$district_id,
  j = 1:C,
  x = 1,
  dims = c(D, C)
)

A <- as.matrix(A_sparse)
district_counts <- rowSums(A)

# stan data list
stan_data <- list(
  C = C,
  D = D,
  y = cluster_df$y,
  n = cluster_df$n,
  district  = cluster_df$district_id,
  A = A,
  district_counts = district_counts
)

###########################################################

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit_iid_null <- stan(
  file    = "iid_null.stan",
  data    = stan_data,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  seed    = 123,
  control = list(
    adapt_delta   = 0.95,    # reduce divergences
    max_treedepth = 12       # allow deeper trees
  )
)

# quick summary of key parameters
print(fit_iid_null,
      pars = c("mu", "sigma_u", "sigma_v"),
      probs = c(0.025, 0.5, 0.975))


fit_iid_null_2 <- stan(
  file    = "iid_null.stan",
  data    = stan_data,
  chains  = 4,
  iter    = 4000,       # doubled from 2000
  warmup  = 2000,       # doubled from 1000
  cores   = 4,
  seed    = 123,
  control = list(
    adapt_delta   = 0.95,
    max_treedepth = 12
  )
)

print(fit_iid_null_2,
      pars  = c("mu", "sigma_u", "sigma_v", "pi_ghana"),
      probs = c(0.025, 0.5, 0.975))

































