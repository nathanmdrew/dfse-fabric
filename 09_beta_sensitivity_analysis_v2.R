# ============================================================
# 09_beta_sensitivity_analysis_v2.R
# Sensitivity analysis: effect of beta factor on model results
# Revisited to address the repeated measures, now handled via GM collapse
# ============================================================

# ---------------------------------
# --- Step 0: Setup             ---
# ---------------------------------

library(tidyverse)
library(EnvStats)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(nlme)
library(rlang)

setwd("C:/Users/vom8/dfse-fabric/")

d <- readRDS(file="data/cleaned_data.RDS")


# -----------------------------------
# --- Step 0: Setup and data prep ---
# -----------------------------------

d <- d |> rename(SampleLocation = `Sample Location`)
  
# Set factor types
d <- d |>
  mutate(
    SampleLocation    = factor(SampleLocation),
    Condition         = factor(Condition),
    Timing_new        = factor(Timing_new),
    ParticipantID     = factor(ParticipantID)
  )

# Collapse Shirt → Sleeve (PI-approved)
d <- d |>
  mutate(
    SampleLocation = fct_recode(SampleLocation, "Sleeve" = "Shirt") |> fct_drop()
  )

# Load MLE betas from script 04
top6_beta <- readRDS("04_output/top6_beta.rds")

# Define beta grid
# Note: MLE betas are per-PAH (~0.11 to 0.21); conventional betas are single values
beta_grid <- tibble(
  beta_label = c("Near-zero (0.01)", "MLE (per-PAH)", "LOD/2 (0.5)", "LOD/√2 (0.707)", "Near-LOD (0.9)"),
  beta_value = list(0.01, top6_beta$beta_mle, 0.5, 1/sqrt(2), 0.9),
  is_mle     = c(FALSE, TRUE, FALSE, FALSE, FALSE)
)

beta_grid

# -----------------------------------------------
# --- Step 1: Impute for each beta    ---
# -----------------------------------------------

impute_with_beta <- function(beta_vals, beta_label, d, top6_beta) {
  if (length(beta_vals) == 1) {
    beta_vals <- rep(beta_vals, nrow(top6_beta))
  }
  stopifnot(length(beta_vals) == nrow(top6_beta))

  for (i in seq_len(nrow(top6_beta))) {
    pc   <- top6_beta$pah_col[i]
    lc   <- top6_beta$lod_col[i]
    nm   <- top6_beta$pah_name[i]
    beta <- beta_vals[i]

    pah_vals <- d[[pc]]
    lod_vals <- d[[lc]]
    is_nd    <- !is.na(pah_vals) & !is.na(lod_vals) & (pah_vals <= lod_vals)

    imp_name <- paste0(nm, "_imp")
    d[[imp_name]] <- pah_vals
    d[[imp_name]][is_nd] <- beta * lod_vals[is_nd]
  }

  imp_cols <- paste0(top6_beta$pah_name, "_imp")

  d |>
    mutate(
      totalPAH_imputed = rowSums(pick(all_of(imp_cols)), na.rm = TRUE),
      log_totalPAH_imp = if_else(totalPAH_imputed > 0, log(totalPAH_imputed), NA_real_),
      beta_label = beta_label
    )
}

# one imputed dataframe per beta setting
imputed_data <- pmap(
  list(beta_grid$beta_value, beta_grid$beta_label),
  \(beta_vals, beta_label) impute_with_beta(beta_vals, beta_label, d, top6_beta)
)
names(imputed_data) <- beta_grid$beta_label

# optional: stacked version for downstream workflows
#imputed_data_all <- list_rbind(imputed_data)


# ==================================================
# Step 1b: GM-collapse repeated measures (by beta set)
# ==================================================

# Geometric mean helper (safe for NA / non-positive values)
gm_safe <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) NA_real_ else exp(mean(log(x)))
}

imputed_data_collapsed <- map(
  imputed_data,
  \(df) {
    imp_cols <- paste0(top6_beta$pah_name, "_imp")

    df |>
      summarize(
        across(all_of(imp_cols), gm_safe),
        .by = c(ParticipantID, Condition, Timing_new, SampleLocation, beta_label)
      ) |>
      mutate(
        totalPAH_imputed = rowSums(pick(all_of(imp_cols)), na.rm = TRUE),
        log_totalPAH_imp = if_else(totalPAH_imputed > 0, log(totalPAH_imputed), NA_real_)
      )
  }
)

names(imputed_data_collapsed) <- names(imputed_data)



# -----------------------------------------------------------
# --- Step 2: Fit GM-collapsed mixed models across betas  ---
# -----------------------------------------------------------

fit_lme_one_beta <- function(df, beta_name = NULL) {
  beta_name <- beta_name %||% dplyr::first(na.omit(df$beta_label))

  df_model <- df |>
    dplyr::filter(!is.na(log_totalPAH_imp), totalPAH_imputed > 0) |>
    dplyr::mutate(
      ParticipantID = as.factor(ParticipantID),
      Condition = as.factor(Condition),
      Timing_new = as.factor(Timing_new),
      SampleLocation = as.factor(SampleLocation)
    )

  tryCatch(
    nlme::lme(
      fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
      random = ~ 1 | ParticipantID,
      weights = nlme::varIdent(form = ~ 1 | SampleLocation),
      data = df_model,
      method = "REML",
      na.action = na.omit,
      control = nlme::lmeControl(opt = "optim")
    ),
    error = \(e) structure(list(error = e$message, beta_label = beta_name), class = "fit_error")
  )
}

extract_lme_fitstats <- function(mod, beta_label) {
  tibble::tibble(
    beta_label = beta_label,
    n_obs = stats::nobs(mod),
    aic = AIC(mod),
    bic = BIC(mod),
    logLik = as.numeric(logLik(mod))
  )
}

extract_lme_fixef <- function(mod, beta_label) {
  tt <- summary(mod)$tTable

  tibble::tibble(
    beta_label = beta_label,
    term = rownames(tt),
    estimate = tt[, "Value"],
    std_error = tt[, "Std.Error"],
    df = tt[, "DF"],
    t_value = tt[, "t-value"],
    p_value = tt[, "p-value"]
  ) |>
    dplyr::mutate(
      t_crit = stats::qt(0.975, df = df),
      conf_low = estimate - t_crit * std_error,
      conf_high = estimate + t_crit * std_error,
      fold_change = exp(estimate),
      pct_change = 100 * (exp(estimate) - 1),
      sig_05 = p_value < 0.05
    ) |>
    dplyr::select(
      beta_label, term, estimate, std_error, df, t_value, p_value, sig_05,
      conf_low, conf_high, fold_change, pct_change
    )
}

# Fit models
models_step2 <- purrr::imap(imputed_data_collapsed, fit_lme_one_beta)

# Keep successful fits for comparison tables
models_step2_ok <- models_step2 |> purrr::keep(~ inherits(.x, "lme"))
models_step2_err <- models_step2 |> purrr::keep(~ inherits(.x, "fit_error"))

step2_fitstats <- purrr::imap_dfr(models_step2_ok, extract_lme_fitstats)
step2_fixef <- purrr::imap_dfr(models_step2_ok, extract_lme_fixef) |>
  dplyr::arrange(term, beta_label)

step2_results <- list(
  models = models_step2,
  models_ok = models_step2_ok,
  model_errors = models_step2_err,
  fitstats = step2_fitstats,
  fixef = step2_fixef
)

step2_results

condition_contrast_results <- purrr::imap_dfr(models_step2_ok, \(mod, label) {
  emmeans::emmeans(mod, "Condition") |>
    emmeans::contrast(method = "pairwise", adjust = "tukey") |>
    as.data.frame() |>
    dplyr::mutate(beta_label = label)
})

beta_order <- beta_grid |>
  dplyr::pull(beta_label)

condition_contrast_results <- condition_contrast_results |>
  dplyr::mutate(
    beta_label = factor(beta_label, levels = beta_order),
    significant = p.value < 0.05
  ) |>
  dplyr::arrange(contrast, beta_label)

condition_contrast_summary <- condition_contrast_results |>
  dplyr::select(beta_label, contrast, estimate, SE, df, t.ratio, p.value, significant)

condition_contrast_summary



# -----------------------------------------------
# --- Step 3: Visualize sensitivity results   ---
# -----------------------------------------------

# Condition contrasts: coefficient plot with 95% CI
p_condition <- condition_contrast_results |>
  dplyr::mutate(
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  ) |>
  ggplot(aes(x = beta_label, y = estimate, color = significant)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "#CC79A7")) +
  facet_wrap(~ contrast, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Sensitivity of Condition Contrasts to Beta Factor",
    subtitle = "Estimates ± 95% CI on log scale (Tukey-adjusted p-values)",
    x = NULL, y = "Estimate (log scale)",
    color = "p < 0.05"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# Timing contrast across beta levels (Doffing vs Donning)
timing_results <- purrr::imap_dfr(models_step2_ok, \(mod, label) {
  emmeans::emmeans(mod, "Timing_new") |>
    emmeans::contrast(method = "pairwise", adjust = "none") |>
    as.data.frame() |>
    dplyr::mutate(beta_label = label)
}) |>
  dplyr::mutate(
    beta_label = factor(beta_label, levels = beta_order),
    significant = p.value < 0.05,
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  )

p_timing <- timing_results |>
  ggplot(aes(x = beta_label, y = estimate, color = significant)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "#CC79A7")) +
  coord_flip() +
  labs(
    title = "Sensitivity of Timing Contrast to Beta Factor",
    subtitle = "Estimate ± 95% CI (no multiplicity adjustment)",
    x = NULL, y = "Estimate (log scale)",
    color = "p < 0.05"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

sensitivity_plot <- p_condition / p_timing +
  patchwork::plot_annotation(title = "Beta Sensitivity Analysis")

sensitivity_plot


# -----------------------------------------------
# --- Step 4: Save artifacts to 09_output      ---
# -----------------------------------------------

out_dir <- "09_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Save data/model/table artifacts as .rds
saveRDS(imputed_data, file.path(out_dir, "imputed_data.rds"))
saveRDS(imputed_data_collapsed, file.path(out_dir, "imputed_data_collapsed.rds"))

saveRDS(models_step2, file.path(out_dir, "models_step2.rds"))
saveRDS(models_step2_ok, file.path(out_dir, "models_step2_ok.rds"))
saveRDS(models_step2_err, file.path(out_dir, "models_step2_err.rds"))

saveRDS(step2_results, file.path(out_dir, "step2_results.rds"))
saveRDS(step2_fitstats, file.path(out_dir, "step2_fitstats.rds"))
saveRDS(step2_fixef, file.path(out_dir, "step2_fixef.rds"))

saveRDS(condition_contrast_results, file.path(out_dir, "condition_contrast_results.rds"))
saveRDS(condition_contrast_summary, file.path(out_dir, "condition_contrast_summary.rds"))
saveRDS(timing_results, file.path(out_dir, "timing_results.rds"))

# Save plot objects as .rds (optional but useful)
saveRDS(p_condition, file.path(out_dir, "p_condition.rds"))
saveRDS(p_timing, file.path(out_dir, "p_timing.rds"))
saveRDS(sensitivity_plot, file.path(out_dir, "sensitivity_plot.rds"))

# Save plots as PNG
ggplot2::ggsave(
  filename = file.path(out_dir, "condition_contrast_sensitivity.png"),
  plot = p_condition,
  width = 10, height = 6, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(out_dir, "timing_contrast_sensitivity.png"),
  plot = p_timing,
  width = 8, height = 5, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(out_dir, "sensitivity_plot_combined.png"),
  plot = sensitivity_plot,
  width = 12, height = 9, dpi = 300
)