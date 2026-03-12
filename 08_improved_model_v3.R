### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-03-12
### Purpose: Further improve the model by dealing with multiple measures per combination
###          of ParticipantID, Sample Location, Condition, and Timing_new. 

# Required packages (install once):
#install.packages("tidyverse")
#install.packages("EnvStats")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("emmeans")
#install.packages("patchwork")
#install.packages("glmmTMB")

library(tidyverse)
library(EnvStats)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(nlme)
library(glmmTMB)

setwd("C:/Users/vom8/dfse-fabric/")

d <- readRDS(file="data/cleaned_data.RDS")

# -----------------------------------
# --- Step 0: Setup and data prep ---
# -----------------------------------


# Set factor types
d <- d |>
  mutate(
    `Sample Location` = factor(`Sample Location`),
    Condition         = factor(Condition),
    Timing_new        = factor(Timing_new),
    ParticipantID     = factor(ParticipantID)
  )

# Collapse Shirt → Sleeve (PI-approved)
d <- d |>
  mutate(
    `Sample Location` = fct_recode(`Sample Location`, "Sleeve" = "Shirt") |> fct_drop()
  )

# Quick verification
str(d[, c("Sample Location", "Condition", "Timing_new", "ParticipantID")])

# Verify new_totalPAH
stopifnot("new_totalPAH" %in% names(d), is.numeric(d$new_totalPAH))

# Check for NAs in modeling variables - should be 0
d |>
  summarise(
    across(
      c(`Sample Location`, Condition, Timing_new, ParticipantID, new_totalPAH),
      \(x) sum(is.na(x))
    )
  )

# Confirm Timing_new levels
levels(d$Timing_new)

d <- d |>
  rename(SampleLocation = `Sample Location`)

# -------------------------------------------------------
# --- Step 1a: Build pair_map and identify top 6 PAHs ---
# -------------------------------------------------------

# Position-dependent - should be okay, but be careful if data structure changes
pah_cols <- 8:22
lod_cols <- 24:38

# Build pair map
pair_map <- tibble(
  pair_id = seq_along(pah_cols),
  pah_col = pah_cols,
  lod_col = lod_cols,
  pah_name = names(d)[pah_cols],
  lod_name = names(d)[lod_cols]
)

# Detection rank: count detects per PAH
pah_detect_rank <- pair_map |>
  mutate(
    n_detect = map2_int(pah_col, lod_col, \(pc, lc) {
      sum(d[[pc]] > d[[lc]], na.rm = TRUE)
    })
  ) |>
  arrange(desc(n_detect))

# Top 6
top6 <- pah_detect_rank |>
  slice_head(n = 6)

top6

# -----------------------------------------------
# --- Step 2: MLE beta-substitution for top 6 ---
# -----------------------------------------------
# Decision point #2: using median LOD per PAH as representative LOD
# Decision point #3: assuming lognormal distribution for each PAH

top6_beta <- top6 |>
  mutate(
    lod_rep = map_dbl(lod_col, \(lc) median(d[[lc]], na.rm = TRUE)),
    mle_fit = pmap(list(pah_col, lod_col), \(pc, lc) {
      pah_vals <- d[[pc]]
      lod_vals <- d[[lc]]
      has_both <- !is.na(pah_vals) & !is.na(lod_vals)
      is_nd    <- has_both & (pah_vals <= lod_vals)

      obs <- if_else(is_nd, lod_vals, pah_vals)[has_both]
      cen <- is_nd[has_both]

      tryCatch(
        elnormCensored(x = obs, censored = cen, method = "mle"),
        error = function(e) {
          warning(paste0("MLE failed: ", e$message))
          NULL
        }
      )
    }),
    mu_hat    = map_dbl(mle_fit, \(f) if (!is.null(f)) f$parameters[["meanlog"]] else NA_real_),
    sigma_hat = map_dbl(mle_fit, \(f) if (!is.null(f)) f$parameters[["sdlog"]] else NA_real_),
    z_hat     = (log(lod_rep) - mu_hat) / sigma_hat,
    beta_mle  = exp(mu_hat + sigma_hat^2 / 2) * pnorm(z_hat - sigma_hat) / pnorm(z_hat) / lod_rep
  ) |>
  select(pair_id, pah_name, pah_col, lod_col, n_detect, lod_rep,
         mu_hat, sigma_hat, z_hat, beta_mle)

top6_beta

# -----------------------------------------------
# --- Step 2b: Impute non-detects for top 6 PAHs
# -----------------------------------------------

# For each top 6 PAH, replace non-detects with beta_mle * LOD
for (i in seq_len(nrow(top6_beta))) {
  pc   <- top6_beta$pah_col[i]
  lc   <- top6_beta$lod_col[i]
  beta <- top6_beta$beta_mle[i]
  nm   <- top6_beta$pah_name[i]

  pah_vals <- d[[pc]]
  lod_vals <- d[[lc]]
  is_nd    <- !is.na(pah_vals) & !is.na(lod_vals) & (pah_vals <= lod_vals)

  # Create new imputed column
  imp_name <- paste0(nm, "_imp")
  d[[imp_name]] <- pah_vals
  d[[imp_name]][is_nd] <- beta * lod_vals[is_nd]
}

# Verify: count remaining zeros per imputed PAH
d |>
  select(ends_with("_imp")) |>
  summarise(across(everything(), \(x) sum(x == 0, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "pah", values_to = "n_zeros")

# ---------------------------------------------------
# --- Step 3: Recompute totalPAH from imputed values
# ---------------------------------------------------

# Sum the 6 imputed PAH columns per row
imp_cols <- names(d) |> str_subset("_imp$")

d <- d |>
  mutate(totalPAH_imputed = rowSums(pick(all_of(imp_cols)), na.rm = TRUE))

# Raw sum of top 6 PAHs only (no imputation)
raw6_cols <- top6_beta$pah_col
d <- d |>
  mutate(totalPAH_raw6 = rowSums(pick(all_of(raw6_cols)), na.rm = TRUE))

# Compare original (all PAH) vs. original (top 6, no imputation) vs. imputed (top 6, MLE)
comparison <- tibble(
  version = c("new_totalPAH (all 15)", "totalPAH_raw6 (top 6, no imputation)", "totalPAH_imputed (top 6, MLE)"),
  n_zero  = c(sum(d$new_totalPAH == 0, na.rm = TRUE),
              sum(d$totalPAH_raw6 == 0, na.rm = TRUE),
              sum(d$totalPAH_imputed == 0, na.rm = TRUE)),
  min     = c(min(d$new_totalPAH, na.rm = TRUE),
              min(d$totalPAH_raw6, na.rm = TRUE),
              min(d$totalPAH_imputed, na.rm = TRUE)),
  median  = c(median(d$new_totalPAH, na.rm = TRUE),
              median(d$totalPAH_raw6, na.rm = TRUE),
              median(d$totalPAH_imputed, na.rm = TRUE)),
  mean    = c(mean(d$new_totalPAH, na.rm = TRUE),
              mean(d$totalPAH_raw6, na.rm = TRUE),
              mean(d$totalPAH_imputed, na.rm = TRUE)),
  max     = c(max(d$new_totalPAH, na.rm = TRUE),
              max(d$totalPAH_raw6, na.rm = TRUE),
              max(d$totalPAH_imputed, na.rm = TRUE))
)

comparison

# No constant needed — all values are positive
d <- d |>
  mutate(log_totalPAH_imp = log(totalPAH_imputed))



# ---------------------------------------------------------
# --- Step 4: Build analysis datasets for 3 strategies  ---
# ---------------------------------------------------------

# 4a) Standardize key columns used in all strategies
sample_loc_col <- if ("SampleLocation" %in% names(d)) "SampleLocation" else "Sample Location"

required_cols <- c("ParticipantID", "Condition", "Timing_new", "totalPAH_imputed", sample_loc_col)
missing_cols <- setdiff(required_cols, names(d))
stopifnot(length(missing_cols) == 0)

d_model_base <- d |>
  transmute(
    ParticipantID,
    Condition,
    Timing_new,
    SampleLocation = .data[[sample_loc_col]],
    totalPAH_imputed
  ) |>
  filter(
    !is.na(ParticipantID),
    !is.na(Condition),
    !is.na(Timing_new),
    !is.na(SampleLocation),
    !is.na(totalPAH_imputed),
    totalPAH_imputed > 0
  ) |>
  mutate(
    ParticipantID = as.factor(ParticipantID),
    Condition = as.factor(Condition),
    Timing_new = as.factor(Timing_new),
    SampleLocation = as.factor(SampleLocation)
  )

# 4b) Strategy A dataset: collapse duplicates with geometric mean
d_gm <- d_model_base |>
  summarize(
    totalPAH_imp = exp(mean(log(totalPAH_imputed))),
    n_rep = n(),
    .by = c(ParticipantID, Condition, SampleLocation, Timing_new)
  ) |>
  mutate(log_totalPAH_imp = log(totalPAH_imp))

# 4c) Strategy B dataset: collapse duplicates with median
d_med <- d_model_base |>
  summarize(
    totalPAH_imp = median(totalPAH_imputed),
    n_rep = n(),
    .by = c(ParticipantID, Condition, SampleLocation, Timing_new)
  ) |>
  mutate(log_totalPAH_imp = log(totalPAH_imp))

# 4d) Strategy C dataset: keep rows, add full-key grouping variable for random effect
d_re <- d_model_base |>
  mutate(
    FullKeyID = interaction(
      ParticipantID, Condition, SampleLocation, Timing_new,
      drop = TRUE, lex.order = TRUE
    ),
    log_totalPAH_imp = log(totalPAH_imputed)
  )

# 4e) Quick size check for A/B/C
step4_sizes <- tibble::tibble(
  dataset = c("d_model_base", "d_gm", "d_med", "d_re"),
  n_rows = c(nrow(d_model_base), nrow(d_gm), nrow(d_med), nrow(d_re))
)

step4_sizes



# -----------------------------------------------------------
# --- Step 5: Fit m5-style mixed models for 3 strategies  ---
# -----------------------------------------------------------

# m5-style fixed effects:
#   log_totalPAH_imp ~ Condition + Timing_new + SampleLocation
# Random intercept for ParticipantID
# Heteroscedastic residual variance by SampleLocation (varIdent)

# 5a) Strategy A: GM-collapsed data
m_gm <- nlme::lme(
  fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
  random = ~ 1 | ParticipantID,
  weights = nlme::varIdent(form = ~ 1 | SampleLocation),
  data = d_gm,
  method = "REML",
  na.action = na.omit,
  control = nlme::lmeControl(opt = "optim")
)

# 5b) Strategy B: Median-collapsed data
m_med <- nlme::lme(
  fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
  random = ~ 1 | ParticipantID,
  weights = nlme::varIdent(form = ~ 1 | SampleLocation),
  data = d_med,
  method = "REML",
  na.action = na.omit,
  control = nlme::lmeControl(opt = "optim")
)

# 5c) Strategy C: Keep all rows + add full-key random effect
m_re <- nlme::lme(
  fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
  random = ~ 1 | ParticipantID / FullKeyID,
  weights = nlme::varIdent(form = ~ 1 | SampleLocation),
  data = d_re,
  method = "REML",
  na.action = na.omit,
  control = nlme::lmeControl(opt = "optim")
)

# Key change: use nlme default optimizer (nlminb) + higher iteration limits
m_re_retry <- nlme::lme(
  fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
  random = ~ 1 | ParticipantID / FullKeyID,
  weights = nlme::varIdent(form = ~ 1 | SampleLocation),
  data = d_re,
  method = "REML",
  na.action = na.omit,
  control = nlme::lmeControl(
    maxIter = 200,
    msMaxIter = 200,
    niterEM = 50,
    tolerance = 1e-6,
    msVerbose = TRUE
  )
)

m_re_retry

# 5d) Compact extraction helpers
extract_fixef <- function(mod, model_name) {
  tt <- summary(mod)$tTable
  tibble::tibble(
    model = model_name,
    term = rownames(tt),
    estimate = tt[, "Value"],
    std_error = tt[, "Std.Error"],
    df = tt[, "DF"],
    t_value = tt[, "t-value"],
    p_value = tt[, "p-value"]
  )
}

extract_fitstats <- function(mod, model_name) {
  tibble::tibble(
    model = model_name,
    n_obs = nobs(mod),
    aic = AIC(mod),
    bic = BIC(mod),
    logLik = as.numeric(logLik(mod))
  )
}

step5_fixef <- dplyr::bind_rows(
  extract_fixef(m_gm, "GM collapse"),
  extract_fixef(m_med, "Median collapse"),
  extract_fixef(m_re_retry, "Random effect for duplicates")
)

step5_fitstats <- dplyr::bind_rows(
  extract_fitstats(m_gm, "GM collapse"),
  extract_fitstats(m_med, "Median collapse"),
  extract_fitstats(m_re_retry, "Random effect for duplicates")
)

list(
  models = list(m_gm = m_gm, m_med = m_med, m_re_retry = m_re_retry),
  step5_fixef = step5_fixef,
  step5_fitstats = step5_fitstats
)



# ------------------------------------------------------------
# --- Step 6: Compare model results across A/B/C strategies ---
# ------------------------------------------------------------

# 6a) Add interval estimates + fold-change interpretation
step6_fixef_comp <- step5_fixef |>
  mutate(
    t_crit = qt(0.975, df = df),
    conf_low = estimate - t_crit * std_error,
    conf_high = estimate + t_crit * std_error,
    fold_change = exp(estimate),
    pct_change = 100 * (exp(estimate) - 1),
    fold_low = exp(conf_low),
    fold_high = exp(conf_high),
    pct_low = 100 * (exp(conf_low) - 1),
    pct_high = 100 * (exp(conf_high) - 1),
    sig_05 = p_value < 0.05
  ) |>
  select(
    model, term, estimate, std_error, df, p_value, sig_05,
    conf_low, conf_high,
    fold_change, fold_low, fold_high,
    pct_change, pct_low, pct_high
  )

# 6b) Focus on key terms (Condition + Timing_new) for direct strategy comparison
step6_key_terms <- step6_fixef_comp |>
  filter(stringr::str_detect(term, "^Condition|^Timing_new")) |>
  arrange(term, model)

# 6c) Optional: location terms separately (often many rows)
step6_location_terms <- step6_fixef_comp |>
  filter(stringr::str_detect(term, "^SampleLocation")) |>
  arrange(term, model)

# 6d) Bring fit stats forward with caution label
# Note: AIC/BIC are less directly comparable when n differs (418 vs 448)
step6_fitstats <- step5_fitstats |>
  mutate(
    n_note = dplyr::case_when(
      n_obs == max(n_obs) ~ "full data (no collapse)",
      TRUE ~ "collapsed data"
    )
  )

# 6e) Compact object for review/saving
step6_results <- list(
  key_terms = step6_key_terms,
  location_terms = step6_location_terms,
  fitstats = step6_fitstats,
  all_terms = step6_fixef_comp
)

step6_results



# ------------------------------------------------------------
# --- Step 7: Primary model decision + sensitivity summary ---
# ------------------------------------------------------------

# 7a) Quantify stability of key inferential terms across strategies
step7_key_stability <- step6_key_terms |>
  mutate(
    direction = dplyr::case_when(
      estimate > 0 ~ "positive",
      estimate < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) |>
  summarize(
    n_models = n(),
    min_estimate = min(estimate, na.rm = TRUE),
    max_estimate = max(estimate, na.rm = TRUE),
    est_range = max_estimate - min_estimate,
    min_fold = min(fold_change, na.rm = TRUE),
    max_fold = max(fold_change, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    max_p = max(p_value, na.rm = TRUE),
    all_significant_0_05 = all(sig_05),
    same_direction = n_distinct(direction) == 1,
    .by = term
  ) |>
  arrange(term)

# 7b) Create model decision table
step7_model_decision <- tibble::tribble(
  ~model, ~role, ~decision_basis,
  "GM collapse", "Primary", "Handles duplicate rows directly; robust to right-skew on imputed positive scale; stable effects vs alternatives.",
  "Median collapse", "Sensitivity", "Checks robustness to non-log-average summarization.",
  "Random effect for duplicates", "Sensitivity", "Retains all rows and models within-fullkey dependence explicitly."
)

# 7c) Final recommendation object
step7_recommendation <- list(
  primary_model = "GM collapse",
  sensitivity_models = c("Median collapse", "Random effect for duplicates"),
  key_term_stability = step7_key_stability,
  fit_stats = step6_fitstats,
  model_decision_table = step7_model_decision
)

step7_recommendation


# ---------------------------------------------------------
# --- Step 8: Visual diagnostics for the 3 mixed models ---
# ---------------------------------------------------------

# Assumes these exist from Step 5:
# m_gm, m_med, m_re_retry

model_list <- list(
  `GM collapse` = m_gm,
  `Median collapse` = m_med,
  `Random effect for duplicates` = m_re_retry
)

make_diag_df <- function(mod, model_name) {
  dat <- nlme::getData(mod)
  resp_var <- all.vars(formula(mod))[1]

  fit0 <- fitted(mod, level = 0)
  rnorm <- residuals(mod, type = "normalized")
  rpear <- residuals(mod, type = "pearson")

  # Align data rows to fitted/residual rows when possible
  rn <- names(fit0)
  if (!is.null(rn) && all(rn %in% rownames(dat))) {
    dat_used <- dat[rn, , drop = FALSE]
  } else {
    dat_used <- dat[seq_along(fit0), , drop = FALSE]
  }

  tibble::tibble(
    model = model_name,
    fitted = as.numeric(fit0),
    resid_norm = as.numeric(rnorm),
    resid_pearson = as.numeric(rpear),
    observed = as.numeric(dat_used[[resp_var]]),
    SampleLocation = as.character(dat_used$SampleLocation)
  ) |>
    dplyr::mutate(
      abs_sqrt_resid = sqrt(abs(resid_norm))
    )
}

diag_df <- purrr::imap_dfr(model_list, make_diag_df)

# 8a) Residuals vs fitted
p_resid_fitted <- ggplot2::ggplot(diag_df, ggplot2::aes(x = fitted, y = resid_norm)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  ggplot2::geom_point(alpha = 0.6, size = 1.6) +
  ggplot2::geom_smooth(se = FALSE, method = "loess", linewidth = 0.7, color = "steelblue4") +
  ggplot2::facet_wrap(~ model, scales = "free_x") +
  ggplot2::labs(
    title = "Residuals vs Fitted",
    x = "Fitted values",
    y = "Normalized residuals"
  ) +
  ggplot2::theme_bw()

# 8b) Normal Q-Q of normalized residuals
p_qq <- ggplot2::ggplot(diag_df, ggplot2::aes(sample = resid_norm)) +
  ggplot2::stat_qq(alpha = 0.6, size = 1.4) +
  ggplot2::stat_qq_line(color = "firebrick", linewidth = 0.7) +
  ggplot2::facet_wrap(~ model, scales = "free") +
  ggplot2::labs(
    title = "Q-Q plot of normalized residuals",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  ggplot2::theme_bw()

# 8c) Scale-location plot
p_scale_location <- ggplot2::ggplot(diag_df, ggplot2::aes(x = fitted, y = abs_sqrt_resid)) +
  ggplot2::geom_point(alpha = 0.6, size = 1.6) +
  ggplot2::geom_smooth(se = FALSE, method = "loess", linewidth = 0.7, color = "darkorange4") +
  ggplot2::facet_wrap(~ model, scales = "free_x") +
  ggplot2::labs(
    title = "Scale-Location",
    x = "Fitted values",
    y = expression(sqrt("|Normalized residuals|"))
  ) +
  ggplot2::theme_bw()

# 8d) Observed vs fitted
p_obs_fitted <- ggplot2::ggplot(diag_df, ggplot2::aes(x = fitted, y = observed)) +
  ggplot2::geom_point(alpha = 0.6, size = 1.6) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  ggplot2::facet_wrap(~ model, scales = "free") +
  ggplot2::labs(
    title = "Observed vs Fitted (log scale response)",
    x = "Fitted values",
    y = "Observed log_totalPAH_imp"
  ) +
  ggplot2::theme_bw()

# 8e) Residual distribution by sample location
p_resid_by_location <- ggplot2::ggplot(
  diag_df,
  ggplot2::aes(x = SampleLocation, y = resid_norm)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  ggplot2::geom_boxplot(outlier.alpha = 0.4) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ model, scales = "free_y") +
  ggplot2::labs(
    title = "Normalized residuals by sample location",
    x = "Sample location",
    y = "Normalized residuals"
  ) +
  ggplot2::theme_bw()

step8_diagnostics <- list(
  diag_df = diag_df,
  p_resid_fitted = p_resid_fitted,
  p_qq = p_qq,
  p_scale_location = p_scale_location,
  p_obs_fitted = p_obs_fitted,
  p_resid_by_location = p_resid_by_location
)

step8_diagnostics



# --------------------------------------------------------
# --- Step 9: Save data frames and plots to 08_output   ---
# --------------------------------------------------------

out_dir <- "08_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 9a) Save all data frame objects currently in environment as .rds
obj_names <- ls(envir = .GlobalEnv)

df_names <- obj_names[
  vapply(obj_names, \(nm) is.data.frame(get(nm, envir = .GlobalEnv)), logical(1))
]

purrr::walk(
  df_names,
  \(nm) {
    saveRDS(
      object = get(nm, envir = .GlobalEnv),
      file = file.path(out_dir, paste0(nm, ".rds"))
    )
  }
)

# 9b) Save diagnostic plots as .png
plot_map <- list(
  p_resid_fitted = p_resid_fitted,
  p_qq = p_qq,
  p_scale_location = p_scale_location,
  p_obs_fitted = p_obs_fitted,
  p_resid_by_location = p_resid_by_location
)

purrr::iwalk(
  plot_map,
  \(p, nm) {
    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(nm, ".png")),
      plot = p,
      width = 10,
      height = 7,
      dpi = 300
    )
  }
)

tibble::tibble(
  saved_dataframes = length(df_names),
  saved_plots = length(plot_map),
  output_dir = out_dir
)



# ---------------------------------------------------------
# --- Step 9b: Save models and key list artifacts (.rds) ---
# ---------------------------------------------------------

# Models
model_artifacts <- list(
  m_gm = m_gm,
  m_med = m_med,
  m_re_retry = m_re_retry
)

purrr::iwalk(
  model_artifacts,
  \(obj, nm) {
    saveRDS(obj, file = file.path(out_dir, paste0(nm, ".rds")))
  }
)

# Key list artifacts from later steps
list_artifacts <- list(
  step6_results = step6_results,
  step7_recommendation = step7_recommendation,
  step8_diagnostics = step8_diagnostics
)

purrr::iwalk(
  list_artifacts,
  \(obj, nm) {
    saveRDS(obj, file = file.path(out_dir, paste0(nm, ".rds")))
  }
)

# Optional manifest of what was saved in this block
step9_saved_objects <- tibble::tibble(
  object = c(names(model_artifacts), names(list_artifacts)),
  type = c(rep("model", length(model_artifacts)), rep("list", length(list_artifacts))),
  path = file.path(out_dir, paste0(c(names(model_artifacts), names(list_artifacts)), ".rds"))
)

step9_saved_objects