###############################################################################
#  PRIOR & POSTERIOR PREDICTIVE CHECKS — BYM2 Fay-Herriot (Ghana 2022)
library(cmdstanr)
library(sf)
library(spdep)
library(INLA)
library(bayesplot)
library(ggplot2)
library(posterior)
library(patchwork)
library(tidyverse)
library(geodata)
library(survey)

set.seed(972)

# PREPARE THE SPATIAL GRAPH
map_sf <- gadm(country = "GHA", level = 2, path = "dhs_gh/") %>%
  st_as_sf() %>%
  st_transform(crs = 4326) %>%
  rename(district_name = NAME_2, district_code = GID_2)

N <- nrow(map_sf)

# Building the neighbours list 
nb <- poly2nb(map_sf)

# Checking for islands (districts with no neighbours)
n_islands <- sum(card(nb) == 0)

# Converting to edge list for Stan (upper triangle only: i < j) 
node1 <- integer(0)
node2 <- integer(0)

for (i in 1:N) {
  neighbors <- nb[[i]]
  if (neighbors[1] != 0) {
    for (j in neighbors) {
      if (i < j) {
        node1 <- c(node1, i)
        node2 <- c(node2, j)
      }
    }
  }
}

N_edges <- length(node1)


# Calculating the BYM2 Scaling Factor
# The scaling factor ensures the ICAR component has unit generalised variance,
# making sigma interpretable as the marginal SD of the combined BYM2 effect
# and phi interpretable as the proportion of variance that is spatial.
# See Riebler et al. (2016), Section 3.

adj_matrix <- nb2mat(nb, style = "B", zero.policy = TRUE)
Q <- diag(rowSums(adj_matrix)) - adj_matrix

# Use INLA's utility to compute the scaled precision matrix
# This applies the sum-to-zero constraint and scales so that
# the geometric mean of the marginal variances equals 1.
Q_scaled <- INLA::inla.scale.model(Q, constr = list(A = matrix(1, 1, N), e = 0))
scaling_factor <- exp(mean(log(diag(Q_scaled))))

# BUILD THE STAN DATA LIST
# Loading your district-level data
df <- read.csv("ghana_data.csv")

# Covariate matrix — ONLY the actual covariates
X_matrix <- as.matrix(df[, c("elevation_std", "rainfall_std", "travel_time_std",
                             "poverty_std", "pop_density_std", "literacy_rate_std")])

# Observed district indices
obs_idx <- which(df$observed == 1)
N_obs   <- length(obs_idx)

df$y_i[is.na(df$y_i)]     <- 0
df$psi_i[is.na(df$psi_i)] <- 1

# --- Assemble the Stan data list ---
stan_data <- list(
  N              = nrow(df),
  K              = ncol(X_matrix),
  y              = df$y_i,
  sigma_e        = df$psi_i,
  X              = X_matrix,
  N_edges        = N_edges,
  node1          = node1,
  node2          = node2,
  scaling_factor = scaling_factor,
  N_obs          = N_obs,
  obs_idx        = obs_idx,
  prior_only     = 1L
)

# RUNING PRIOR PREDICTIVE CHECK (prior_only = 1)
# Compile the Stan model
mod <- cmdstan_model("prior_predictive_checks.stan")

fit_prior <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 500,
  seed            = 2026
)

# Check that sampling worked
fit_prior$cmdstan_diagnose()


#

# Extract draws
draws <- fit_prior$draws(format = "matrix")

# --- 4a. Prior predictive prevalences ---
prev_draws <- fit_prior$draws("prev", format = "matrix")   # [n_draws x N]
y_rep_draws <- fit_prior$draws("y_rep", format = "matrix")  # [n_draws x N]

cat("\n=== PRIOR PREDICTIVE DIAGNOSTIC SUMMARY ===\n")
cat("% district prevalences in [1%, 60%]:",
    round(100 * mean(prev_draws > 0.01 & prev_draws < 0.60), 1), "%\n")
cat("% district prevalences < 1%:",
    round(100 * mean(prev_draws < 0.01), 1), "%\n")
cat("% district prevalences > 60%:",
    round(100 * mean(prev_draws > 0.60), 1), "%\n")
cat("% district prevalences > 90%:",
    round(100 * mean(prev_draws > 0.90), 1), "%\n")


# --- PLOT 1: Distribution of implied district prevalences ---
p1 <- tibble(prev = as.vector(prev_draws)) |>
  ggplot(aes(x = prev)) +
  geom_histogram(bins = 100, fill = "#2171b5", colour = "white",
                 linewidth = 0.1) +
  geom_vline(xintercept = c(0.01, 0.60), linetype = "dashed",
             colour = "red", linewidth = 0.5) +
  annotate("text", x = 0.30, y = Inf, vjust = 2,
           label = "Plausible range for\nunder-5 malaria in W. Africa",
           colour = "red", size = 3) +
  labs(
    title = "Prior Predictive: Implied District Prevalences",
    subtitle = "All draws × all districts — does this look plausible?",
    x = "Prevalence", y = "Count"
  ) +
  theme_minimal(base_size = 12)


# --- PLOT 2: 50 random prior predictive prevalence densities ---
n_overlay <- 50
sim_idx <- sample(1:nrow(prev_draws), n_overlay)

p2 <- tibble(
  sim  = rep(sim_idx, each = N),
  prev = as.vector(t(prev_draws[sim_idx, ]))
) |>
  ggplot(aes(x = prev, group = factor(sim))) +
  geom_density(alpha = 0.2, fill = "#2171b5", colour = "#08519c",
               linewidth = 0.3) +
  labs(
    title = paste0(n_overlay, " Prior Predictive Draws"),
    subtitle = "Each curve = one simulated Ghana under the prior",
    x = "Prevalence", y = "Density"
  ) +
  theme_minimal(base_size = 12)


# --- PLOT 3: Prior predictive for y_rep vs observed y ---
p3 <- ppc_dens_overlay(
  y    = stan_data$y,
  yrep = y_rep_draws[1:min(100, nrow(y_rep_draws)), ]
) +
  labs(
    title = "Prior Predictive Check: y_rep vs Observed y",
    subtitle = "Dark line = observed direct estimates; light lines = simulated from priors"
  )


# --- PLOT 4: Spatial heterogeneity (IQR of prevalences per draw) ---
iqr_per_draw <- apply(prev_draws, 1, function(x) diff(quantile(x, c(0.25, 0.75))))

p4 <- tibble(iqr = iqr_per_draw) |>
  ggplot(aes(x = iqr)) +
  geom_histogram(bins = 50, fill = "#d94801", colour = "white",
                 linewidth = 0.1) +
  labs(
    title = "Prior Predictive: Spatial Heterogeneity",
    subtitle = "IQR of district prevalences within each simulation",
    x = "IQR", y = "Count"
  ) +
  theme_minimal(base_size = 12)



combined <- (p1 + p2) / ( p4) +
  plot_annotation(
    title = "Prior Predictive Checks — BYM2 Fay-Herriot, Ghana 2022",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(combined)
ggsave("prior_predictive_checks.png", combined, width = 14, height = 10, dpi = 300)


# --- PLOT 5: Prior distributions of key parameters ---
sigma_draws <- fit_prior$draws("sigma", format = "matrix")
phi_draws   <- fit_prior$draws("phi", format = "matrix")
beta0_draws <- fit_prior$draws("beta0", format = "matrix")

p5a <- mcmc_hist(data.frame(sigma = as.vector(sigma_draws)), binwidth = 0.02) +
  labs(title = "Prior on σ (BYM2 total SD)", subtitle = "PC: P(σ > 1) = 0.01")

p5b <- mcmc_hist(data.frame(phi = as.vector(phi_draws)), binwidth = 0.02) +
  labs(title = "Prior on φ (spatial fraction)", subtitle = "Beta(0.5, 0.5)")

p5c <- mcmc_hist(data.frame(beta0 = as.vector(beta0_draws)), binwidth = 0.05) +
  labs(title = "Prior on β₀ (intercept)", subtitle = "N(-2, 1)")

(p5a + p5b + p5c) +
  plot_annotation(title = "Prior Distributions of Key Parameters")

ggsave("prior_distributions.png", width = 14, height = 4, dpi = 300)



