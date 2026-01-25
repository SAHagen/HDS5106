library(tidyverse)
library(sf)         
library(haven)      
library(survey)  
library(tmap)
library(srvyr)
library(terra)
library(exactextractr)

# reading in the data

who_boundaries <- read_sf("glob_polygons/GLOBAL_ADM2.shp") # for country boundaries (polygons)
dhs_clusters <- read_sf("dhs_gh/GHGE8AFL.shp") # dhs coordinates  
dhs_children <- read_dta("dhs_gh/GHKR8CFL.dta") # Kr file for ACT treatment and Fever variables
dhs_hm <- read_dta("dhs_gh/GHPR8CFL.dta") # PR file for nets (ITNs)and malaria test results

#Filtering for Ghana polygons and re projecting them  into the WGS 84 coordinate reference system
ghana_districts <- who_boundaries %>%
  filter(ADM0_NAME == "GHANA") %>%
  st_transform(crs = 4326)

# cluster locations from DHS GPS file to join to the polygons above
gps <- dhs_clusters %>% 
  rename(cluster_id = DHSCLUST)



gps_joined <- st_join(gps, ghana_districts, join = st_nearest_feature)

# preparing the KR file for ACt treatment and fever variables
kr <- dhs_children %>% 
  select(
    # v001 is Cluster number (Primary Sampling Unit). Used for merging and survey design.
    cluster_id = v001,
    # v002 is the household number used for merging within the cluster
    hh_id = v002,
    #b16 is the line number of the child in the household,  
    line_id = b16,
    #h22 had fever in the last 2 week or not
    fever = h22,
    #Child took ACT for fever in last 2 weeks.
    took_act = h32z
  ) %>% 
  mutate(
    act_binary = case_when(
      took_act == 1 ~ 1,
      took_act == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

#household member PR data for nets (ITNs) and malaria test results 
pr <- dhs_hm %>% 
  select(cluster_id = hv001,
         hh_id = hv002,
         line_id = hvidx,
         weight_raw = hv005,
         stratum = hv023,
         age_months = hc1,
         malaria_results = hml35,
         net_use = hml12,
         wealth = hv270,
         urban_rural = hv025) %>% 
  filter(age_months < 60) %>% 
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

#joining the PR and KR file
my_data_dhs <- kr %>%
  left_join(pr, by = c("cluster_id", "hh_id", "line_id"))

#joing the PR, KR files to the Polygons data
my_df <- my_data_dhs %>%
  left_join(gps_joined, by = c("cluster_id"))

#creating the survey design that correctly reflects the DHS sampling scheme.

# This accounts for: 
# 1) Clustering (individuals within the same PSU/cluster are correlated), 
# 2) Stratification (sampling differed by strata), and 
# 3) Sampling weights (to make estimates representative of the national population,
#   where urban observations represent ~1:5000 people and rural observations ~1:500).
# This ensures all estimates (means, prevalence, CIs) calcuated from this part are unbiased and have
# correct standard errors.

survey_design <- my_df %>%
    filter(!is.na(weight)) %>% # no weights are calculated for missing oberservations, so we filter them out
  as_survey_design(
    ids = cluster_id,   
    strata = stratum,   
    weights = weight    
  )

# calculating the district level estimates
district_estimates <- survey_design %>%
  group_by(ADM2_NAME) %>% 
    summarise(malaria_prevalence = survey_mean(outcome_binary, vartype = "ci", na.rm = TRUE),
              act_coverage = survey_mean(act_binary, vartype = "ci", na.rm = TRUE),
              itn_coverage = survey_mean(itn_binary, vartype = "ci", na.rm = TRUE))

# joining the polygons with estimates fro ploting.
map_data_polygons <- ghana_districts %>%
  left_join(district_estimates , by = "ADM2_NAME") %>%
  mutate(
    # Creating the hover label on the polygon object
    hover_label = paste0(ADM1_NAME, ": ", round(malaria_prevalence * 100, 1), "%")
  )


tmap_mode("view")

tm_shape(map_data_polygons) +
  tm_polygons(
    col = "prevalence",           
    style = "pretty",
    palette = "-RdBu", 
    title = "Malaria Prevelance",
    # hover settings 
    id = "hover_label",           
    popup.vars = c(               
      "State" = "ADM1_NAME", 
      "Rate" = "prevalence",
      "Lower CI" = "prevalence_low",
      "Upper CI" = "prevalence_upp"
    )
  )

act_data_polygons <- ghana_districts %>%
  left_join(district_estimates , by = "ADM2_NAME") %>%
  mutate(
    # Creating the hover label on the polygon object
    hover_label = paste0(ADM1_NAME, ": ", round(act_coverage * 100, 1), "%")
  )

act_cov <- tm_shape(act_data_polygons) +
  tm_polygons(
    col = "prevalence",           
    style = "pretty",
    palette = "Greens", 
    title = "ACT upatake",
    # hover settings 
    id = "hover_label",           
    popup.vars = c(               
      "State" = "ADM1_NAME", 
      "Rate" = "prevalence",
      "Lower CI" = "prevalence_low",
      "Upper CI" = "prevalence_upp"
    )
  )


itn_data_polygons <- ghana_districts %>%
  left_join(district_estimates , by = "ADM2_NAME") %>%
  mutate(
    # Creating the hover label on the polygon object
    hover_label = paste0(ADM1_NAME, ": ", round(itn_coverage * 100, 1), "%")
  )

itn_cov <- tm_shape(itn_data_polygons) +
  tm_polygons(
    col = "prevalence",           
    style = "pretty",
    palette = "Purples", 
    title = "ACT upatake",
    # hover settings 
    id = "hover_label",           
    popup.vars = c(               
      "State" = "ADM1_NAME", 
      "Rate" = "prevalence",
      "Lower CI" = "prevalence_low",
      "Upper CI" = "prevalence_upp"
    )
  )


act_cov
itn_cov
malaria_prev

