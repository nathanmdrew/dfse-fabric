### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-03-13
### Purpose: Simulate beta-substitution under a favorable censoring scenario (~30% ND)
###          using MLE-based beta estimation (Ganser & Hewett approach)

library(tidyverse)
library(EnvStats)

setwd("C:/Users/vom8/dfse-fabric/")

# ============================================================
# Beta-Substitution Simulation: Single Draw
# Favorable scenario: ~30% non-detect, continuous LOD
# ============================================================

set.seed(716)
n <- 500

# --- 1. True lognormal PAH concentrations ---
# Parameters chosen to loosely resemble a detectable PAH
mu_log    <- 1.0   # meanlog
sigma_log <- 0.8   # sdlog
true_vals <- rlnorm(n, meanlog = mu_log, sdlog = sigma_log)

# --- 2. Continuous LOD distribution (lognormal, lower/narrower than PAH) ---
# lod_meanlog tuned so that ~30% of true values fall below their paired LOD
lod_meanlog <- 0.5
lod_sdlog   <- 0.4
lod_vals    <- rlnorm(n, meanlog = lod_meanlog, sdlog = lod_sdlog)

# --- 3. Classify: detect vs. non-detect ---
is_detect <- true_vals > lod_vals
is_nd     <- !is_detect

cat("Empirical censoring rate:", round(100 * mean(is_nd), 1), "%\n")

# --- 4. Build observed vector as a lab would record it ---
# Detected: true value is recorded.
# Non-detect: only the LOD is known (the censoring point).
obs_vals <- if_else(is_nd, lod_vals, true_vals)

# --- 5. Estimate beta via MLE of censored lognormal (Ganser & Hewett) ---
# elnormCensored() fits a lognormal to the full censored dataset via MLE.
# censored = TRUE flags left-censored observations (value recorded at LOD).
# Beta is derived from the truncated lognormal expectation below lod_ref:
#   beta = E[X | X <= lod_ref] / lod_ref
fit <- elnormCensored(
  x        = obs_vals,
  censored = is_nd,
  method   = "mle"
)

mu_hat    <- fit$parameters[["meanlog"]]
sigma_hat <- fit$parameters[["sdlog"]]
lod_ref   <- median(lod_vals)   # representative LOD for beta derivation

cat("MLE mu_hat:", round(mu_hat, 3), "| sigma_hat:", round(sigma_hat, 3), "\n")

# Truncated lognormal beta derivation
z        <- (log(lod_ref) - mu_hat) / sigma_hat
beta_hat <- exp(mu_hat + sigma_hat^2 / 2) * pnorm(z - sigma_hat) / pnorm(z) / lod_ref

cat("Estimated beta (MLE):", round(beta_hat, 3), "\n")

# --- 6. Apply beta-substitution to non-detects ---
imputed_vals <- if_else(is_nd, beta_hat * lod_vals, true_vals)

# --- 7. Assemble results tibble ---
sim_dat <- tibble(
  true_val    = true_vals,
  lod_val     = lod_vals,
  obs_val     = obs_vals,
  is_detect   = is_detect,
  is_nd       = is_nd,
  imputed_val = imputed_vals
)

# ============================================================
# Evaluation Metric 1: GM, GSD, and P95 recovery
# ============================================================
gm_fn  <- function(x) exp(mean(log(x[x > 0])))
gsd_fn <- function(x) exp(sd(log(x[x > 0])))

metrics <- tribble(
  ~metric, ~true,                             ~imputed,
  "GM",    gm_fn(sim_dat$true_val),           gm_fn(sim_dat$imputed_val),
  "GSD",   gsd_fn(sim_dat$true_val),          gsd_fn(sim_dat$imputed_val),
  "P95",   quantile(sim_dat$true_val, 0.95),  quantile(sim_dat$imputed_val, 0.95)
) |>
  mutate(pct_error = 100 * (imputed - true) / true)

metrics

# ============================================================
# Evaluation Metric 2: KDE overlay - True vs. Beta-Substituted
# ============================================================
plot_sim <- sim_dat |>
  select(true_val, imputed_val) |>
  pivot_longer(everything(), names_to = "distribution", values_to = "value") |>
  mutate(
    distribution = recode(distribution,
      "true_val"    = "True (full)",
      "imputed_val" = "Beta-Substituted"
    )
  )

p_kde <- ggplot(plot_sim, aes(x = value, fill = distribution, color = distribution)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values  = c("True (full)"      = "#E69F00",
                                "Beta-Substituted" = "#0072B2")) +
  scale_color_manual(values = c("True (full)"      = "#E69F00",
                                "Beta-Substituted" = "#0072B2")) +
  labs(
    title    = "Simulation: True vs. Beta-Substituted Distribution (MLE beta)",
    subtitle = paste0(
      round(100 * mean(is_nd), 1), "% censoring | n = ", n,
      " | beta = ",     round(beta_hat, 3),
      " | GM error = ", round(metrics$pct_error[metrics$metric == "GM"], 1), "%"
    ),
    x = "Concentration", y = "Density",
    fill = NULL, color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

p_kde

# ============================================================
# Evaluation Metric 3: Scatter - True vs. Imputed (ND only)
# ============================================================
p_scatter <- sim_dat |>
  filter(is_nd) |>
  ggplot(aes(x = true_val, y = imputed_val)) +
  geom_point(alpha = 0.4, color = "#0072B2") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title    = "Beta-Substituted vs. True Value (Non-detects only)",
    subtitle = "Dashed line = perfect recovery | MLE-based beta",
    x = "True Concentration",
    y = "Beta-Substituted Concentration"
  ) +
  theme_minimal(base_size = 13)

p_scatter
