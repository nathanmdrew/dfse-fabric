  true_val    = true_vals,
  lod_val     = lod_vals,
  obs_val     = obs_vals,
  is_detect   = is_detect, By convention, censored observations are recorded at the censoring limit (LOD).
# Beta is derived from the truncated lognormal expectation below the median LOD:
# beta = E[X | X <= LOD_ref] / LOD_ref
fit <- elnormCensored(
  x        = obs_vals,  # detected values + LOD for non-detects
  censored = is_nd,     # TRUE = left-censored at LOD
  method   = "mle"
)

# EnvStats returns parameters on the log scale named "mean" and "sd"
mu_hat    <- fit$parameters[["mean"]]
sigma_hat <- fit$parameters[["sd"]]
lod_ref   <- median(lod_vals)  # representative LOD for beta derivation

cat("MLE mu_hat:", round(mu_hat, 3), "| sigma_hat:", round(sigma_hat, 3), "\n")

# Beta from truncated lognormal: E[X | X <= LOD_ref] / LOD_ref
z        <- (log(lod_ref) - mu_hat) / sigma_hat
beta_hat <- exp(mu_hat + sigma_hat^2 / 2) * pnorm(z - sigma_hat) / pnorm(z) / lod_ref

cat("Estimated beta (MLE):", round(beta_hat, 3), "\n")nlog
sigma_log <- 0.8   # sdlog
true_vals <- rlnorm(n, meanlog = mu_log, sdlog = sigma_log)

# --- 2. Continuous LOD distribution (lognormal, lower/narrower than PAH) ---
# Tune lod_meanlog so ~30% of true values fall below their paired LOD
lod_meanlog <- 0.5
lod_sdlog   <- 0.4
lod_vals    <- rlnorm(n, meanlog = lod_meanlog, sdlog = lod_sdlog)

# --- 3. Classify: detect vs. non-detect ---
is_detect <- true_vals > lod_vals
is_nd     <- !is_detect

cat("Empirical censoring rate:", round(100 * mean(is_nd), 1), "%\n")

# --- 4. Estimate beta via MLE of censored lognormal (Ganser & Hewett) ---
# elnormCensored() fits a lognormal to the full censored dataset (detects +
# non-detects). Beta is then derived from the truncated lognormal expectation:
# E[X | X <= LOD] / LOD, evaluated at the median LOD as a representative
# threshold for the continuous LOD distribution.
fit <- elnormCensored(
  x              = true_vals,
  censored       = is_nd,       # TRUE = non-detect (left-censored)
  censoring.side = "left",
  method         = "mle"
)

mu_hat    <- fit$parameters[["meanlog"]]
sigma_hat <- fit$parameters[["sdlog"]]
lod_ref   <- median(lod_vals)   # representative LOD for beta derivation

# Beta from truncated lognormal: E[X | X <= LOD] / LOD
z        <- (log(lod_ref) - mu_hat) / sigma_hat
beta_hat <- exp(mu_hat + sigma_hat^2 / 2) * pnorm(z - sigma_hat) / pnorm(z) / lod_ref

cat("MLE mu_hat:", round(mu_hat, 3), "| sigma_hat:", round(sigma_hat, 3), "\n")
cat("Estimated beta (MLE):", round(beta_hat, 3), "\n")

# --- 5. Apply beta-substitution to non-detects ---
imputed_vals <- if_else(is_nd, beta_hat * lod_vals, true_vals)

# --- 6. Assemble results tibble ---
sim_dat <- tibble(
  true_val    = true_vals,
  lod_val     = lod_vals,
  is_detect   = is_detect,
  is_nd       = is_nd,
  imputed_val = imputed_vals
)

# ============================================================
# Evaluation Metric 1: GM, GSD, and P95 recovery
# ============================================================
gm  <- function(x) exp(mean(log(x[x > 0])))
gsd <- function(x) exp(sd(log(x[x > 0])))

metrics <- tribble(
  ~metric, ~true,                             ~imputed,
  "GM",    gm(sim_dat$true_val),              gm(sim_dat$imputed_val),
  "GSD",   gsd(sim_dat$true_val),             gsd(sim_dat$imputed_val),
  "P95",   quantile(sim_dat$true_val, 0.95),  quantile(sim_dat$imputed_val, 0.95)
) |>
  mutate(pct_error = 100 * (imputed - true) / true)

metrics

# ============================================================
# Evaluation Metric 2: KDE overlay - True vs. β-Substituted
# ============================================================
plot_sim <- sim_dat |>
  select(true_val, imputed_val) |>
  pivot_longer(everything(), names_to = "distribution", values_to = "value") |>
  mutate(
    distribution = recode(distribution,
      "true_val"    = "True (full)",
      "imputed_val" = "β-Substituted"
    )
  )

p_kde <- ggplot(plot_sim, aes(x = value, fill = distribution, color = distribution)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values  = c("True (full)"   = "#E69F00",
                                "β-Substituted" = "#0072B2")) +
  scale_color_manual(values = c("True (full)"   = "#E69F00",
                                "β-Substituted" = "#0072B2")) +
  labs(
    title    = "Simulation: True vs. β-Substituted Distribution (MLE beta)",
    subtitle = paste0(
      "~30% censoring | n = ", n,
      " | β̂ = ",       round(beta_hat, 3),
      " | GM error = ", round(filter(metrics, metric == "GM")$pct_error, 1), "%"
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
    title    = "β-Substituted vs. True Value (Non-detects only)",
    subtitle = "Dashed line = perfect recovery | MLE-based beta",
    x = "True Concentration",
    y = "β-Substituted Concentration"
  ) +
  theme_minimal(base_size = 13)

p_scatter
