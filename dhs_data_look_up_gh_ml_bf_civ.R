# data preparation for multi year dhs sae

# Question?? Are there similar variables across the different countries and years
# checking for Mali, CIV and Burkina Faso and Ghana. If yes which are the tests recorded for each countries
# This helps us decide whether to have the country on board for modelling or drop it
library(haven)
library(sf)
library(tidyverse)

survey_index <- tribble(
  ~country, ~iso,  ~year, ~pr_path,                              ~gps_path,
  "Ghana",  "GHA", 2022,  "dhs_gh_22/GHPR8CDT/GHPR8CFL.DTA",    "dhs_gh_22/GHGE8AFL/GHGE8AFL.shp",
  "Ghana",  "GHA", 2014,  "dhs_gh_14/GHPR72DT/GHPR72FL.DTA",    "dhs_gh_14/GHGE71FL/GHGE71FL.shp",
  "Burkina","BFA", 2021,  "dhs_bf_21/BFPR81DT/BFPR81FL.DTA",    "dhs_bf_21/BFGE81FL/BFGE81FL.shp",
  "Burkina","BFA", 2010,  "dhs_bf_10/BFPR62DT/BFPR62FL.DTA",    "dhs_bf_10/BFGE61FL/BFGE61FL.shp",
  "Mali",   "MLI", 2023,  "dhs_ml_23/MLPR8ADT/MLPR8AFL.dta",    "dhs_ml_23/MLGE8AFL/MLGE8AFL.shp",
  "Mali",   "MLI", 2018,  "dhs_ml_18/MLPR7ADT/MLPR7AFL.DTA",    "dhs_ml_18/MLGE7AFL/MLGE7AFL.shp",
  "Mali",   "MLI", 2012,  "dhs_ml_12/MLPR6ADT/MLPR6AFL.DTA",    "dhs_ml_12/MLGE6BFL/MLGE6BFL.shp",
  "CIV",    "CIV", 2021,  "dhs_ci_21/CIPR81DT/CIPR81FL.DTA",    "dhs_ci_21/CIGE81FL/CIGE81FL.shp",
  "CIV",    "CIV", 2012,  "dhs_ci_12/CIPR62DT/CIPR62FL.DTA",    "dhs_ci_12/CIGE61FL/CIGE61FL.shp"
)


check_survey <- function(pr_path, gps_path, country, year) {
  cat("\n=========================================\n")
  cat(country, year, "\n")
  cat("=========================================\n")
  
  df <- read_dta(pr_path, n_max = 200)
  
  # Key variables
  vars_needed <- c("hv001", "hv005", "hv023", "hc1", "hml32", "hml35")
  present     <- vars_needed %in% tolower(names(df))
  names(present) <- vars_needed
  
  cat("Variable check:\n")
  for (v in vars_needed) {
    cat(" ", v, ":", ifelse(present[v], "PRESENT", "MISSING"), "\n")
  }
  
  # RDT coding if present
  if (present["hml32"]) {
    col <- df[[which(tolower(names(df)) == "hml32")]]
    cat("hml32 values:", paste(sort(unique(as.integer(col))), collapse = ", "), "\n")
  }
  
  # GPS cluster variable
  gps <- read_sf(gps_path)
  cat("GPS cluster var (DHSCLUST):", "DHSCLUST" %in% names(gps), "\n")
  cat("GPS rows (clusters):", nrow(gps), "\n")
}

# Runing for all surveys
purrr::pwalk(
  survey_index,
  ~check_survey(..4, ..5, ..1, ..3)
)



check_mali <- function(pr_path, year) {
  cat("\n=== Mali", year, "===\n")
  df <- read_dta(pr_path)
  
  # Full table of hml32 values
  cat("hml32 full table:\n")
  print(table(df$hml32, useNA = "always"))
  
  # Full table of hml35 values
  cat("hml35 full table:\n")
  print(table(df$hml35, useNA = "always"))
  
  # How many children under 5?
  under5 <- df[!is.na(df$hc1) & df$hc1 < 60, ]
  cat("Children under 5:", nrow(under5), "\n")
  
  # Of those, how many have valid hml32?
  cat("Under 5 with hml32 = 0 or 1:",
      sum(under5$hml32 %in% c(0, 1), na.rm = TRUE), "\n")
  
  # Of those, how many have valid hml35?
  cat("Under 5 with hml35 = 0 or 1:",
      sum(under5$hml35 %in% c(0, 1), na.rm = TRUE), "\n")
}

check_mali("dhs_ml_18/MLPR7ADT/MLPR7AFL.DTA", 2018)
check_mali("dhs_ml_23/MLPR8ADT/MLPR8AFL.dta", 2023)

check_civ <- function(pr_path, year) {
  cat("\n=== CIV", year, "===\n")
  df <- read_dta(pr_path)
  
  cat("hml32 full table:\n")
  print(table(df$hml32, useNA = "always"))
  
  # Check labels
  cat("\nhml32 value labels:\n")
  print(attr(df$hml32, "labels"))
  
  # Children under 5
  under5 <- df[!is.na(df$hc1) & df$hc1 < 60, ]
  cat("\nUnder 5 — hml32 table:\n")
  print(table(under5$hml32, useNA = "always"))
}

check_civ("dhs_ci_21/CIPR81DT/CIPR81FL.DTA", 2021)


# Mali has no Malaria results for their latest dhs round and the two available have different test for both years
# therefor cannot be used because of comparability purposes.