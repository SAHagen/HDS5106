library(tidyverse)
library(sf)          # For spatial handling
library(haven)       # To read DHS .dta files 
library(survey)     # For survey weighted analysis 
library(tmap)
library(srvyr)
library(terra)
library(exactextractr)

# boundaries from the WHO dataset
who_boundaries <- read_sf("GLOBAL_ADM2.shp")
nig_states <- who_boundaries %>%
  filter(ADM0_NAME == "NIGERIA") %>%
  st_transform(crs = 4326)

# loading DHS Data Points
dhs_clusters <- read_sf("NGGE8AFL.shp") %>%
  rename(cluster_id = DHSCLUST)


joined_data <- st_join(dhs_clusters, nig_states, join = st_nearest_feature)

# Calculating wieghted avarages
dhs_children <- read_dta("NGKR8AFL.dta") %>%
   select(v001, v002, v005, v021, v022, h9) %>%  #v001 clsuter number, v002 household number, v005 sample wieght, v021 sampling unit, v022 sample strate, h9 measles
   mutate(
     # Creating the weight variable scaled by 1000000
     wt = v005 / 1000000,
     
     # Creating binary variable: 1 if vaccinated/has measles, 0 if not
     has_measles_vax = case_when(
       # GROUP 1: CONFIRMED YES
       # 1 = Date on card, 2 = Mother said yes, 3 = Marked on card
       h9 %in% c(1, 2, 3) ~ 1,  
       
       # GROUP 2: CONFIRMED NO
       h9 == 0 ~ 0,             
       
       # GROUP 3: HANDLING "DON'T KNOW" (8) , will treat these NA vales
       h9 == 8 ~ NA_real_,
       
       TRUE ~ NA_real_
     )
   )

my_data <- dhs_children %>%
  left_join(joined_data, by = c("v001" = "cluster_id"))


survey_design <- my_data %>%
  as_survey_design(
    ids = v021,      # Cluster ID
    strata = v022,   # Stratification
    weights = wt     # Calculated Weight
  )


state_estimates <- survey_design %>%
  group_by(ADM1_NAME) %>% 
    summarise(prevalence = survey_mean(has_measles_vax, vartype = "ci", na.rm = TRUE))
  


###Ploting#######
tmap_mode("view")
map_data_interactive <- map_data %>%
  mutate(
    hover_label = paste0(ADM2_NAME, ": ", round(prevalence * 100, 1), "%")
  )


tm_shape(map_data_interactive) +
  tm_polygons(
    col = "prevalence",           
    style = "pretty",
    palette = "RdYlBu", 
    title = "Vaccination Rate",
    
    # INTERACTIVITY SETTINGS:
    id = "hover_label",           
    popup.vars = c(               
      "State" = "ADM1_NAME", 
      "Rate" = "prevalence",
      "Lower CI" = "prevalence_low",
      "Upper CI" = "prevalence_upp"
    )
  )




