# ============================================================
# 12_main_analysis.R
# Implement main decision points
# Author: Nathan M. Drew
# Date: 2026-03-16
# ============================================================

# ---------------------------------
# --- Step 0: Setup             ---
# ---------------------------------

library(tidyverse)
library(EnvStats)
library(emmeans)
library(patchwork)
library(nlme) #for heteroscedastic mixed model

setwd("C:/Users/vom8/dfse-fabric/")

d <- readRDS(file="data/cleaned_data.RDS")

d <- d |> rename(SampleLocation = `Sample Location`)

# Collapse Shirt → Sleeve (PI-approved)
d <- d |>
  mutate(
    SampleLocation = fct_recode(SampleLocation, "Sleeve" = "Shirt") |> fct_drop()
  )


# ---------------------------------
# --- Step 1: Define Top 4 PAHs  ---
# ---------------------------------

# Goal: Identify and freeze the primary response definition (Top 4 PAHs by global detect rate).
# Logic: Use raw measured data + matching LOD columns, compute detect rates once, rank, and select top 4.

# 1a) Define PAH and LOD column positions from cleaned dataset structure
pah_cols <- 8:22
lod_cols <- 24:38

stopifnot(length(pah_cols) == 15, length(lod_cols) == 15)

pah_lookup_15 <- tibble::tibble(
  pah_name = names(d)[pah_cols],
  pah_col = names(d)[pah_cols],
  lod_col = names(d)[lod_cols]
)

# 1b) Compute global detect summary (raw measured values; no imputation here)
detect_summary <- purrr::pmap_dfr(
  pah_lookup_15,
  \(pah_name, pah_col, lod_col) {
    x <- d[[pah_col]]
    lod <- d[[lod_col]]

    n_nonmissing <- sum(!is.na(x) & !is.na(lod))
    n_detect <- sum(!is.na(x) & !is.na(lod) & x > lod)

    tibble::tibble(
      pah_name = pah_name,
      pah_col = pah_col,
      lod_col = lod_col,
      n_nonmissing = n_nonmissing,
      n_detect = n_detect,
      pct_detect = dplyr::if_else(n_nonmissing > 0, 100 * n_detect / n_nonmissing, NA_real_)
    )
  }
) |>
  dplyr::arrange(dplyr::desc(pct_detect)) |>
  dplyr::mutate(rank_detect = dplyr::row_number())

# 1c) Freeze Top 4 set for primary analysis
top4_map <- detect_summary |>
  dplyr::slice_head(n = 4) |>
  dplyr::select(rank_detect, pah_name, pah_col, lod_col, n_detect, pct_detect)

top4_pah_names <- top4_map$pah_name
top4_pah_cols <- top4_map$pah_col
top4_lod_cols <- top4_map$lod_col

# 1d) Tiny QC checks for reproducibility
step1_qc <- list(
  n_total_pahs = nrow(detect_summary),
  n_top4 = length(top4_pah_cols),
  top4_names = top4_pah_names,
  top4_is_sorted = identical(top4_map$rank_detect, 1:4),
  paired_cols_same_length = length(top4_pah_cols) == length(top4_lod_cols)
)

stopifnot(
  step1_qc$n_total_pahs == 15,
  step1_qc$n_top4 == 4,
  step1_qc$top4_is_sorted,
  step1_qc$paired_cols_same_length
)

step1_artifacts <- list(
  pah_lookup_15 = pah_lookup_15,
  detect_summary = detect_summary,
  top4_map = top4_map,
  top4_pah_names = top4_pah_names,
  top4_pah_cols = top4_pah_cols,
  top4_lod_cols = top4_lod_cols,
  step1_qc = step1_qc
)

step1_artifacts



# ---------------------------------
# --- Step 2: MLE imputation     ---
# ---------------------------------

# Goal: Impute non-detects for the Top 4 PAHs using the same PAH-level MLE beta approach as script 11.
# Logic: Fit censored lognormal per PAH (EnvStats::elnormCensored), derive one beta per PAH, impute ND as beta * LOD_i.

# 2a) Helper: convert fitted lognormal params to PAH-level beta at a reference LOD
calc_beta_from_fit <- function(meanlog, sdlog, lod_ref) {
  z1 <- (log(lod_ref) - meanlog) / sdlog
  z2 <- (log(lod_ref) - meanlog - sdlog^2) / sdlog
  ex_trunc <- exp(meanlog + 0.5 * sdlog^2) * stats::pnorm(z2) / stats::pnorm(z1)
  ex_trunc / lod_ref
}

# 2b) Helper: fit one PAH + impute censored values
impute_one_pah_envstats <- function(x, lod) {
  ok <- !is.na(x) & !is.na(lod)
  x_ok <- x[ok]
  lod_ok <- lod[ok]

  censored <- x_ok <= lod_ok
  n_total <- length(x_ok)
  n_cens <- sum(censored)
  n_detect <- n_total - n_cens

  x_imp <- x
  beta <- NA_real_
  meanlog <- NA_real_
  sdlog <- NA_real_
  fit_method <- NA_character_
  fit_ok <- FALSE

  # Edge cases
  if (n_total == 0) {
    return(list(
      x_imp = x_imp, beta = beta, meanlog = meanlog, sdlog = sdlog,
      fit_method = "no_complete_pairs", fit_ok = FALSE,
      n_total = n_total, n_detect = n_detect, n_cens = n_cens
    ))
  }

  if (n_detect == 0) {
    x_imp[ok] <- 0.5 * lod_ok
    return(list(
      x_imp = x_imp, beta = 0.5, meanlog = meanlog, sdlog = sdlog,
      fit_method = "all_censored_fallback_half_lod", fit_ok = FALSE,
      n_total = n_total, n_detect = n_detect, n_cens = n_cens
    ))
  }

  if (n_cens == 0) {
    return(list(
      x_imp = x_imp, beta = NA_real_, meanlog = NA_real_, sdlog = NA_real_,
      fit_method = "no_censoring", fit_ok = TRUE,
      n_total = n_total, n_detect = n_detect, n_cens = n_cens
    ))
  }

  # Fit censored lognormal (use LOD as censoring limit for censored rows)
  x_fit <- x_ok
  x_fit[censored] <- lod_ok[censored]

  fit <- try(
    EnvStats::elnormCensored(
      x = x_fit,
      censored = censored,
      method = "mle"
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    x_ok_imp <- x_ok
    x_ok_imp[censored] <- 0.5 * lod_ok[censored]
    x_imp[ok] <- x_ok_imp

    return(list(
      x_imp = x_imp, beta = 0.5, meanlog = meanlog, sdlog = sdlog,
      fit_method = "fit_failed_fallback_half_lod", fit_ok = FALSE,
      n_total = n_total, n_detect = n_detect, n_cens = n_cens
    ))
  }

  meanlog <- unname(fit$parameters[["meanlog"]])
  sdlog <- unname(fit$parameters[["sdlog"]])

  # Script-11-consistent reference LOD: median of all non-missing LOD values
  lod_ref <- stats::median(lod_ok, na.rm = TRUE)
  if (!is.finite(lod_ref) || lod_ref <= 0) {
    lod_ref <- stats::median(lod_ok[lod_ok > 0], na.rm = TRUE)
  }

  beta <- calc_beta_from_fit(meanlog = meanlog, sdlog = sdlog, lod_ref = lod_ref)

  # Impute censored observations only
  x_ok_imp <- x_ok
  x_ok_imp[censored] <- beta * lod_ok[censored]
  x_imp[ok] <- x_ok_imp

  list(
    x_imp = x_imp, beta = beta, meanlog = meanlog, sdlog = sdlog,
    fit_method = "elnormCensored_mle_beta_times_lod", fit_ok = TRUE,
    n_total = n_total, n_detect = n_detect, n_cens = n_cens
  )
}

# 2c) Apply imputation across Top 4 PAHs
d_imp <- d
imp_log_list <- vector("list", length(top4_pah_cols))
top4_imp_cols <- paste0(top4_pah_cols, "_imp")

for (i in seq_along(top4_pah_cols)) {
  pc <- top4_pah_cols[i]
  lc <- top4_lod_cols[i]
  ic <- top4_imp_cols[i]

  res <- impute_one_pah_envstats(d_imp[[pc]], d_imp[[lc]])
  d_imp[[ic]] <- res$x_imp

  imp_log_list[[i]] <- tibble::tibble(
    pah_name = top4_map$pah_name[i],
    pah_col = pc,
    lod_col = lc,
    imp_col = ic,
    beta_mle = res$beta,
    meanlog = res$meanlog,
    sdlog = res$sdlog,
    fit_method = res$fit_method,
    fit_ok = res$fit_ok,
    n_total = res$n_total,
    n_detect = res$n_detect,
    n_cens = res$n_cens
  )
}

imputation_log <- dplyr::bind_rows(imp_log_list)

# 2d) Build imputed total PAH response (row-level, pre-collapse)
d_imp <- d_imp |>
  dplyr::mutate(
    n_nonmissing_imp = rowSums(!is.na(dplyr::pick(dplyr::all_of(top4_imp_cols)))),
    totalPAH_imputed = rowSums(dplyr::pick(dplyr::all_of(top4_imp_cols)), na.rm = TRUE),
    totalPAH_imputed = dplyr::if_else(n_nonmissing_imp == 0L, NA_real_, totalPAH_imputed)
  )

# 2e) Tiny QC for Step 2
step2_qc <- list(
  n_top4_pahs = length(top4_pah_cols),
  n_imp_cols = length(top4_imp_cols),
  all_fit_ok = all(imputation_log$fit_ok),
  any_fallback = any(!imputation_log$fit_ok),
  min_beta = min(imputation_log$beta_mle, na.rm = TRUE),
  max_beta = max(imputation_log$beta_mle, na.rm = TRUE)
)

step2_artifacts <- list(
  d_imp = d_imp,
  top4_imp_cols = top4_imp_cols,
  imputation_log = imputation_log,
  step2_qc = step2_qc
)

step2_artifacts



# ---------------------------------
# --- Step 3: GM-collapse        ---
# ---------------------------------

# Goal: Collapse repeated measurements within ParticipantID × Condition × Timing_new × SampleLocation.
# Logic: Use geometric mean (GM) on imputed Top-4 PAH columns, then rebuild totalPAH and log-response at collapsed level.

# 3a) Helper: safe geometric mean for positive values
gm_safe <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) NA_real_ else exp(mean(log(x)))
}

# 3b) Collapse to one row per modeling unit and GM-collapse each imputed PAH
d_gm <- d_imp |>
  dplyr::summarize(
    dplyr::across(dplyr::all_of(top4_imp_cols), gm_safe),
    .by = c(ParticipantID, Condition, Timing_new, SampleLocation)
  )

# 3c) Recompute response variables at collapsed level
d_gm <- d_gm |>
  dplyr::mutate(
    n_nonmissing_imp = rowSums(!is.na(dplyr::pick(dplyr::all_of(top4_imp_cols)))),
    totalPAH_imputed = rowSums(dplyr::pick(dplyr::all_of(top4_imp_cols)), na.rm = TRUE),
    totalPAH_imputed = dplyr::if_else(n_nonmissing_imp == 0L, NA_real_, totalPAH_imputed),
    log_totalPAH_imp = dplyr::if_else(totalPAH_imputed > 0, log(totalPAH_imputed), NA_real_)
  )

# 3d) Tiny QC checks
step3_qc <- list(
  n_rows_pre = nrow(d_imp),
  n_rows_post = nrow(d_gm),
  n_collapsed = nrow(d_imp) - nrow(d_gm),
  n_response_nonmissing = sum(!is.na(d_gm$log_totalPAH_imp)),
  expected_post_rows = d_imp |>
    dplyr::distinct(ParticipantID, Condition, Timing_new, SampleLocation) |>
    nrow()
)

stopifnot(step3_qc$n_rows_post == step3_qc$expected_post_rows)

step3_artifacts <- list(
  d_gm = d_gm,
  step3_qc = step3_qc
)

step3_artifacts



# ---------------------------------
# --- Step 4: Fit primary LME    ---
# ---------------------------------

# Goal: Fit the primary heteroscedastic mixed model on the GM-collapsed dataset.
# Logic: Use additive fixed effects for Condition, Timing_new, and SampleLocation;
#        random intercept for ParticipantID; allow residual variance to differ by SampleLocation.

# 4a) Prepare modeling dataframe
d_model <- d_gm |>
  dplyr::filter(!is.na(log_totalPAH_imp)) |>
  dplyr::mutate(
    ParticipantID = as.factor(ParticipantID),
    Condition = as.factor(Condition),
    Timing_new = as.factor(Timing_new),
    SampleLocation = as.factor(SampleLocation)
  ) |>
  droplevels()

# 4b) Fit primary heteroscedastic model (no fallback in primary script)
warn_step4 <- character(0)

m_primary <- withCallingHandlers(
  nlme::lme(
    fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
    random = ~ 1 | ParticipantID,
    weights = nlme::varIdent(form = ~ 1 | SampleLocation),
    data = d_model,
    method = "REML",
    na.action = na.omit,
    control = nlme::lmeControl(opt = "optim")
  ),
  warning = \(w) {
    warn_step4 <<- c(warn_step4, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

# 4c) Tiny QC for fit object + sample size
step4_qc <- tibble::tibble(
  n_rows_model = nrow(d_model),
  n_participants = dplyr::n_distinct(d_model$ParticipantID),
  n_condition = dplyr::n_distinct(d_model$Condition),
  n_timing = dplyr::n_distinct(d_model$Timing_new),
  n_location = dplyr::n_distinct(d_model$SampleLocation),
  aic = AIC(m_primary),
  bic = BIC(m_primary),
  logLik = as.numeric(logLik(m_primary)),
  n_warnings = length(warn_step4)
)

step4_artifacts <- list(
  d_model = d_model,
  m_primary = m_primary,
  warnings = warn_step4,
  step4_qc = step4_qc
)

step4_artifacts



# ----------------------------------------------
# --- Step 4b: Model diagnostics + plots      ---
# ----------------------------------------------

# Goal: Check basic model assumptions and identify obvious outliers/influence patterns.
# Logic: Use fitted values and normalized residuals from nlme model, then create standard diagnostic plots.

# 4b-1) Build diagnostic dataframe
diag_df <- d_model |>
  dplyr::mutate(
    fitted = stats::fitted(m_primary),
    resid_raw = stats::residuals(m_primary, type = "response"),
    resid_norm = stats::residuals(m_primary, type = "normalized"),
    abs_resid_norm = abs(resid_norm)
  )

# 4b-2) Simple diagnostic summaries
diag_summary <- tibble::tibble(
  n_obs = nrow(diag_df),
  resid_norm_mean = mean(diag_df$resid_norm, na.rm = TRUE),
  resid_norm_sd = stats::sd(diag_df$resid_norm, na.rm = TRUE),
  n_abs_resid_gt2 = sum(diag_df$abs_resid_norm > 2, na.rm = TRUE),
  n_abs_resid_gt3 = sum(diag_df$abs_resid_norm > 3, na.rm = TRUE),
  pct_abs_resid_gt3 = 100 * mean(diag_df$abs_resid_norm > 3, na.rm = TRUE)
)

# 4b-3) Plot: residuals vs fitted
p_diag_resid_fitted <- ggplot2::ggplot(diag_df, ggplot2::aes(x = fitted, y = resid_norm)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_point(alpha = 0.7, size = 1.8) +
  ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.7, color = "steelblue") +
  ggplot2::labs(
    title = "Residuals vs Fitted",
    x = "Fitted values",
    y = "Normalized residuals"
  ) +
  ggplot2::theme_minimal(base_size = 11)

# 4b-4) Plot: normal Q-Q of normalized residuals
qq_obj <- stats::qqnorm(diag_df$resid_norm, plot.it = FALSE)
qq_df <- tibble::tibble(theoretical = qq_obj$x, sample = qq_obj$y)

p_diag_qq <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
  ggplot2::geom_point(alpha = 0.7, size = 1.6) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  ggplot2::labs(
    title = "Q-Q Plot (Normalized Residuals)",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  ggplot2::theme_minimal(base_size = 11)

# 4b-5) Plot: residual spread by SampleLocation (checks variance pattern)
p_diag_resid_by_loc <- ggplot2::ggplot(diag_df, ggplot2::aes(x = SampleLocation, y = abs_resid_norm)) +
  ggplot2::geom_boxplot(outlier.alpha = 0.4) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "|Normalized residuals| by SampleLocation",
    x = NULL,
    y = "|Normalized residual|"
  ) +
  ggplot2::theme_minimal(base_size = 10)

# 4b-6) Optional combined diagnostics panel
p_diag_combined <- p_diag_resid_fitted / p_diag_qq / p_diag_resid_by_loc

# Add to Step 4 artifacts
step4_artifacts <- c(
  step4_artifacts,
  list(
    diag_df = diag_df,
    diag_summary = diag_summary,
    p_diag_resid_fitted = p_diag_resid_fitted,
    p_diag_qq = p_diag_qq,
    p_diag_resid_by_loc = p_diag_resid_by_loc,
    p_diag_combined = p_diag_combined
  )
)

diag_summary



# ---------------------------------
# --- Step 5: Primary inference  ---
# ---------------------------------

# Goal: Extract model-based estimates for manuscript reporting.
# Logic: Summarize fixed effects, fit stats, and pairwise contrasts for Condition, Timing_new, and SampleLocation.

# 5a) Model fit statistics
fitstats_step5 <- tibble::tibble(
  n_obs = stats::nobs(m_primary),
  aic = AIC(m_primary),
  bic = BIC(m_primary),
  logLik = as.numeric(logLik(m_primary))
)

# 5b) Fixed-effects table (log scale + exponentiated interpretation)
tt <- summary(m_primary)$tTable

fixef_step5 <- tibble::tibble(
  term = rownames(tt),
  estimate = tt[, "Value"],
  std_error = tt[, "Std.Error"],
  df = tt[, "DF"],
  t_value = tt[, "t-value"],
  p_value = tt[, "p-value"]
) |>
  dplyr::mutate(
    fold_change = exp(estimate),
    pct_change = 100 * (exp(estimate) - 1)
  )

# 5c) Condition contrasts (primary target; Tukey adjustment)
condition_contrast_results <- emmeans::emmeans(m_primary, "Condition") |>
  emmeans::contrast(method = "pairwise", adjust = "tukey") |>
  as.data.frame() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    significant = p.value < 0.05
  )

# 5d) Timing contrast (single pair; no multiplicity adjustment)
timing_results <- emmeans::emmeans(m_primary, "Timing_new") |>
  emmeans::contrast(method = "pairwise", adjust = "none") |>
  as.data.frame() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    significant = p.value < 0.05
  )

# 5e) Sampling-location contrasts (PI request; many comparisons, Tukey-adjusted)
location_contrast_results <- emmeans::emmeans(m_primary, "SampleLocation") |>
  emmeans::contrast(method = "pairwise", adjust = "tukey") |>
  as.data.frame() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    significant = p.value < 0.05
  )

# 5f) Bundle Step 5 artifacts
step5_artifacts <- list(
  fitstats = fitstats_step5,
  fixef = fixef_step5,
  condition_contrast_results = condition_contrast_results,
  timing_results = timing_results,
  location_contrast_results = location_contrast_results
)

step5_artifacts



# ---------------------------------
# --- Step 6: Primary figures    ---
# ---------------------------------

# Goal: Create manuscript-ready contrast plots for Condition, Timing_new, and SampleLocation.
# Logic: Use estimated contrasts from Step 5, add 95% CI, and draw a zero-reference line.

# 6a) Add 95% CI columns to each contrast table
condition_plot_df <- condition_contrast_results |>
  dplyr::mutate(
    ci_low = estimate - stats::qt(0.975, df = df) * SE,
    ci_high = estimate + stats::qt(0.975, df = df) * SE
  )

timing_plot_df <- timing_results |>
  dplyr::mutate(
    ci_low = estimate - stats::qt(0.975, df = df) * SE,
    ci_high = estimate + stats::qt(0.975, df = df) * SE
  )

location_plot_df <- location_contrast_results |>
  dplyr::mutate(
    ci_low = estimate - stats::qt(0.975, df = df) * SE,
    ci_high = estimate + stats::qt(0.975, df = df) * SE
  )

# 6b) Condition contrast plot (primary)
p_condition_step6 <- condition_plot_df |>
  ggplot2::ggplot(ggplot2::aes(x = contrast, y = estimate, color = significant)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_low, ymax = ci_high), size = 0.6) +
  ggplot2::scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey50")) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Condition contrasts (Top 4 PAHs, MLE-imputed, GM-collapsed)",
    x = NULL,
    y = "Estimate (log scale)",
    color = "p < 0.05"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "top")

# 6c) Timing contrast plot
p_timing_step6 <- timing_plot_df |>
  ggplot2::ggplot(ggplot2::aes(x = contrast, y = estimate, color = significant)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_low, ymax = ci_high), size = 0.8) +
  ggplot2::scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "grey50")) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Timing contrast",
    x = NULL,
    y = "Estimate (log scale)",
    color = "p < 0.05"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "top")

# 6d) Location contrasts plot (secondary/PI-requested)
p_location_step6 <- location_plot_df |>
  ggplot2::ggplot(ggplot2::aes(x = estimate, y = forcats::fct_reorder(contrast, estimate), color = significant)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = ci_low, xmax = ci_high), size = 0.35) +
  ggplot2::scale_color_manual(values = c("TRUE" = "#CC79A7", "FALSE" = "grey60")) +
  ggplot2::labs(
    title = "SampleLocation pairwise contrasts (Tukey-adjusted)",
    x = "Estimate (log scale)",
    y = NULL,
    color = "p < 0.05"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "top")

# 6e) Combined panel for quick review
p_step6_combined <- p_condition_step6 / p_timing_step6

step6_artifacts <- list(
  condition_plot_df = condition_plot_df,
  timing_plot_df = timing_plot_df,
  location_plot_df = location_plot_df,
  p_condition_step6 = p_condition_step6,
  p_timing_step6 = p_timing_step6,
  p_location_step6 = p_location_step6,
  p_step6_combined = p_step6_combined
)

step6_artifacts



# ---------------------------------
# --- Step 7: Save artifacts     ---
# ---------------------------------

# Goal: Save all primary-analysis artifacts for reproducibility and manuscript drafting.
# Logic: Write analysis objects as .rds and figures as .png into 12_output.

output_dir_12 <- "12_output"
dir.create(output_dir_12, showWarnings = FALSE, recursive = TRUE)

# 7a) Save step artifact bundles (.rds)
saveRDS(step1_artifacts, file.path(output_dir_12, "step1_artifacts.rds"))
saveRDS(step2_artifacts, file.path(output_dir_12, "step2_artifacts.rds"))
saveRDS(step3_artifacts, file.path(output_dir_12, "step3_artifacts.rds"))
saveRDS(step4_artifacts, file.path(output_dir_12, "step4_artifacts.rds"))
saveRDS(step5_artifacts, file.path(output_dir_12, "step5_artifacts.rds"))
saveRDS(step6_artifacts, file.path(output_dir_12, "step6_artifacts.rds"))

# 7b) Save key standalone objects (.rds)
saveRDS(detect_summary, file.path(output_dir_12, "detect_summary.rds"))
saveRDS(top4_map, file.path(output_dir_12, "top4_map.rds"))
saveRDS(imputation_log, file.path(output_dir_12, "imputation_log.rds"))
saveRDS(m_primary, file.path(output_dir_12, "m_primary.rds"))

saveRDS(fitstats_step5, file.path(output_dir_12, "fitstats_step5.rds"))
saveRDS(fixef_step5, file.path(output_dir_12, "fixef_step5.rds"))
saveRDS(condition_contrast_results, file.path(output_dir_12, "condition_contrast_results.rds"))
saveRDS(timing_results, file.path(output_dir_12, "timing_results.rds"))
saveRDS(location_contrast_results, file.path(output_dir_12, "location_contrast_results.rds"))

# 7c) Save figures (.png)
ggplot2::ggsave(
  filename = file.path(output_dir_12, "condition_contrasts.png"),
  plot = p_condition_step6,
  width = 9, height = 5, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_12, "timing_contrast.png"),
  plot = p_timing_step6,
  width = 8, height = 4.5, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_12, "location_contrasts.png"),
  plot = p_location_step6,
  width = 10, height = 12, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_12, "condition_timing_combined.png"),
  plot = p_step6_combined,
  width = 10, height = 8, dpi = 300
)

# 7d) Save manifest for quick audit
step7_manifest <- tibble::tibble(
  file = list.files(output_dir_12, full.names = FALSE)
) |>
  dplyr::arrange(file)

saveRDS(step7_manifest, file.path(output_dir_12, "step7_manifest.rds"))

step7_manifest

# 7e) Save diagnostic objects + plots from Step 4b
saveRDS(diag_df, file.path(output_dir_12, "diag_df.rds"))
saveRDS(diag_summary, file.path(output_dir_12, "diag_summary.rds"))
saveRDS(p_diag_resid_fitted, file.path(output_dir_12, "p_diag_resid_fitted.rds"))
saveRDS(p_diag_qq, file.path(output_dir_12, "p_diag_qq.rds"))
saveRDS(p_diag_resid_by_loc, file.path(output_dir_12, "p_diag_resid_by_loc.rds"))
saveRDS(p_diag_combined, file.path(output_dir_12, "p_diag_combined.rds"))

ggplot2::ggsave(
  filename = file.path(output_dir_12, "diag_resid_fitted.png"),
  plot = p_diag_resid_fitted,
  width = 8, height = 4.5, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_12, "diag_qq.png"),
  plot = p_diag_qq,
  width = 6, height = 5, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_12, "diag_resid_by_location.png"),
  plot = p_diag_resid_by_loc,
  width = 8, height = 7, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_12, "diag_combined.png"),
  plot = p_diag_combined,
  width = 10, height = 12, dpi = 300
)



# ---------------------------------
# --- Step 8 (optional): Handoff ---
# ---------------------------------

# Goal: Save a lightweight handoff record for coauthors and future reproducibility checks.
# Logic: Store key analysis decisions + session metadata as portable artifacts.

# 8a) Key decision log
analysis_decisions <- tibble::tibble(
  decision_id = c(
    "response_definition",
    "non_detect_imputation",
    "repeat_measure_handling",
    "primary_model",
    "condition_contrasts",
    "timing_contrast",
    "location_contrasts"
  ),
  decision = c(
    "Top 4 PAHs (global detect-rate ranking from raw measured data)",
    "PAH-level MLE beta via EnvStats::elnormCensored; ND imputed as beta * LOD_i",
    "GM collapse within ParticipantID x Condition x Timing_new x SampleLocation",
    "nlme::lme with random intercept (ParticipantID) and varIdent(~1|SampleLocation), REML",
    "emmeans pairwise contrasts with Tukey adjustment",
    "emmeans pairwise contrast with adjust='none' (single contrast)",
    "emmeans pairwise contrasts with Tukey adjustment (secondary/PI-requested)"
  )
)

# 8b) Session metadata snapshot
session_meta <- list(
  analysis_date = Sys.time(),
  r_version = R.version.string,
  package_versions = tibble::tibble(
    package = c("tidyverse", "EnvStats", "nlme", "emmeans", "patchwork"),
    version = c(
      as.character(packageVersion("tidyverse")),
      as.character(packageVersion("EnvStats")),
      as.character(packageVersion("nlme")),
      as.character(packageVersion("emmeans")),
      as.character(packageVersion("patchwork"))
    )
  ),
  script = "12_main_analysis.R"
)

# 8c) Save optional handoff artifacts
saveRDS(analysis_decisions, file.path(output_dir_12, "analysis_decisions.rds"))
saveRDS(session_meta, file.path(output_dir_12, "session_meta.rds"))

# Optional plain-text session info file
writeLines(capture.output(sessionInfo()), con = file.path(output_dir_12, "sessionInfo.txt"))

step8_handoff <- list(
  analysis_decisions = analysis_decisions,
  session_meta = session_meta
)

step8_handoff