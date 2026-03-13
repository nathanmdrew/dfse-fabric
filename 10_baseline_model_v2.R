# ============================================================
# 10_baseline_model_v2.R
# Revisit the baseline models (no imputation)
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

d <- d |> rename(SampleLocation = `Sample Location`)

# Collapse Shirt → Sleeve (PI-approved)
d <- d |>
  mutate(
    SampleLocation = fct_recode(SampleLocation, "Sleeve" = "Shirt") |> fct_drop()
  )

# Reuse-ready paths
artifact_dir_09 <- "09_output"
output_dir_10 <- "10_output"
dir.create(output_dir_10, showWarnings = FALSE, recursive = TRUE)

# Optional reusable artifacts from script 09 (load only if present)
reuse_files <- c(
  "step2_results.rds",
  "condition_contrast_results.rds",
  "timing_results.rds"
)

reuse_paths <- file.path(artifact_dir_09, reuse_files)
reuse_exists <- file.exists(reuse_paths)
names(reuse_exists) <- reuse_files

reuse_09 <- list(
  step2_results = if (reuse_exists["step2_results.rds"]) readRDS(file.path(artifact_dir_09, "step2_results.rds")) else NULL,
  condition_contrast_results = if (reuse_exists["condition_contrast_results.rds"]) readRDS(file.path(artifact_dir_09, "condition_contrast_results.rds")) else NULL,
  timing_results = if (reuse_exists["timing_results.rds"]) readRDS(file.path(artifact_dir_09, "timing_results.rds")) else NULL
)

# Required analysis columns for baseline modeling
required_cols <- c("ParticipantID", "Condition", "Timing_new", "SampleLocation")
missing_cols <- setdiff(required_cols, names(d))
stopifnot(length(missing_cols) == 0)

# Standardize modeling fields (kept consistent across all runs)
d <- d |>
  mutate(
    ParticipantID = as.factor(ParticipantID),
    Condition = as.factor(Condition),
    Timing_new = as.factor(Timing_new),
    SampleLocation = as.factor(SampleLocation)
  )

# Step 0 artifact object
step0_setup <- list(
  output_dir_10 = output_dir_10,
  artifact_dir_09 = artifact_dir_09,
  reuse_exists = reuse_exists,
  n_rows = nrow(d),
  n_cols = ncol(d)
)

step0_setup


# ---------------------------------
# --- Step 1: Define PAH sets    ---
# ---------------------------------

# Column positions provided
pah_cols <- 8:22
lod_cols <- 24:38

stopifnot(length(pah_cols) == 15, length(lod_cols) == 15)

# Build lookup from existing data structure
pah_lookup_15 <- tibble(
  pah_name = names(d)[pah_cols],
  pah_col = names(d)[pah_cols],
  lod_col = names(d)[lod_cols]
)

# Global detect-rate summary (fixed once)
detect_summary <- purrr::pmap_dfr(
  pah_lookup_15,
  \(pah_name, pah_col, lod_col) {
    x <- d[[pah_col]]
    lod <- d[[lod_col]]

    n_nonmissing <- sum(!is.na(x) & !is.na(lod))
    n_detect <- sum(!is.na(x) & !is.na(lod) & x > lod)

    tibble(
      pah_name = pah_name,
      pah_col = pah_col,
      lod_col = lod_col,
      n_nonmissing = n_nonmissing,
      n_detect = n_detect,
      pct_detect = if_else(n_nonmissing > 0, 100 * n_detect / n_nonmissing, NA_real_)
    )
  }
) |>
  arrange(desc(pct_detect))

# Freeze PAH sets globally
pah_sets <- list(
  top1 = detect_summary |> slice_head(n = 1) |> pull(pah_name),
  top4 = detect_summary |> slice_head(n = 4) |> pull(pah_name),
  top6 = detect_summary |> slice_head(n = 6) |> pull(pah_name),
  top15 = detect_summary |> slice_head(n = 15) |> pull(pah_name)
)

pah_sets_map <- purrr::imap(
  pah_sets,
  \(set_names, set_id) {
    detect_summary |>
      filter(pah_name %in% set_names) |>
      mutate(set_id = set_id) |>
      select(set_id, pah_name, pah_col, lod_col, pct_detect)
  }
)

step1_artifacts <- list(
  pah_lookup_15 = pah_lookup_15,
  detect_summary = detect_summary,
  pah_sets = pah_sets,
  pah_sets_map = pah_sets_map
)

saveRDS(step1_artifacts, file.path(output_dir_10, "step1_artifacts.rds"))

# ---------------------------------
# --- Step 1a: QC checks         ---
# ---------------------------------

# Hard checks (stop if violated)
stopifnot(length(pah_cols) == 15)
stopifnot(length(lod_cols) == 15)
stopifnot(all(pah_cols %in% seq_along(d)))
stopifnot(all(lod_cols %in% seq_along(d)))
stopifnot(length(intersect(pah_cols, lod_cols)) == 0)
stopifnot(!anyDuplicated(pah_cols))
stopifnot(!anyDuplicated(lod_cols))

pah_is_numeric <- purrr::map_lgl(names(d)[pah_cols], \(nm) is.numeric(d[[nm]]))
lod_is_numeric <- purrr::map_lgl(names(d)[lod_cols], \(nm) is.numeric(d[[nm]]))
stopifnot(all(pah_is_numeric), all(lod_is_numeric))

# Pair-level QC table (review before pipeline run)
pair_qc <- tibble(
  idx = seq_along(pah_cols),
  pah_col = names(d)[pah_cols],
  lod_col = names(d)[lod_cols]
) |>
  mutate(
    n_complete_pair = purrr::map2_int(pah_col, lod_col, \(pc, lc) sum(!is.na(d[[pc]]) & !is.na(d[[lc]]))),
    n_detect = purrr::map2_int(pah_col, lod_col, \(pc, lc) sum(!is.na(d[[pc]]) & !is.na(d[[lc]]) & d[[pc]] > d[[lc]])),
    pct_detect = if_else(n_complete_pair > 0, 100 * n_detect / n_complete_pair, NA_real_),
    pct_non_detect = if_else(n_complete_pair > 0, 100 - pct_detect, NA_real_)
  ) |>
  mutate(
    qc_flag = case_when(
      n_complete_pair == 0 ~ "NO_COMPLETE_PAIR_DATA",
      n_detect == 0 ~ "ZERO_DETECTS",
      TRUE ~ "OK"
    )
  )

# Detect-summary consistency check
detect_summary_qc <- detect_summary |>
  mutate(
    pct_detect_in_range = dplyr::between(pct_detect, 0, 100),
    rank = row_number()
  )

stopifnot(all(detect_summary_qc$pct_detect_in_range, na.rm = TRUE))

# Top-set membership check
pah_sets_check <- tibble(
  set = c("top1", "top4", "top6", "top15"),
  n = c(length(pah_sets$top1), length(pah_sets$top4), length(pah_sets$top6), length(pah_sets$top15))
)
stopifnot(identical(pah_sets_check$n, c(1L, 4L, 6L, 15L)))

step1_qc <- list(
  pair_qc = pair_qc,
  detect_summary_qc = detect_summary_qc,
  pah_sets_check = pah_sets_check,
  flagged_pairs = pair_qc |> filter(qc_flag != "OK")
)

saveRDS(step1_qc, file.path(output_dir_10, "step1_qc.rds"))

step1_qc


# ---------------------------------
# --- Step 1b: Analysis run grid ---
# ---------------------------------

analysis_grid <- tidyr::crossing(
  pah_set = c("top1", "top4", "top6", "top15"),
  response_mode = c("all_data_plus_c", "detects_only")
) |>
  dplyr::mutate(
    run_id = dplyr::row_number(),
    run_label = paste0("run_", run_id, "_", pah_set, "_", response_mode),
    n_pahs = dplyr::case_when(
      pah_set == "top1" ~ 1L,
      pah_set == "top4" ~ 4L,
      pah_set == "top6" ~ 6L,
      pah_set == "top15" ~ 15L
    ),
    pah_names = purrr::map(pah_set, \(s) pah_sets[[s]])
  ) |>
  dplyr::relocate(run_id, run_label, pah_set, n_pahs, response_mode, pah_names)

saveRDS(analysis_grid, file.path(output_dir_10, "analysis_grid.rds"))

analysis_grid

# ---------------------------------
# --- Step 1c: Build run datasets ---
# ---------------------------------

build_run_dataset <- function(run_id, run_label, pah_set, response_mode, pah_names, d, detect_summary) {
  set_map <- detect_summary |>
    dplyr::filter(pah_name %in% pah_names) |>
    dplyr::select(pah_name, pah_col, lod_col)

  pah_value_cols <- set_map$pah_col
  lod_value_cols <- set_map$lod_col

  gm_nonzero_or_zero <- function(x) {
    x_obs <- x[!is.na(x)]
    if (length(x_obs) == 0) return(NA_real_)
    x_pos <- x_obs[x_obs > 0]
    if (length(x_pos) == 0) return(0)
    exp(mean(log(x_pos)))
  }

  # row-level detect flag (Option A basis)
  detect_mat <- purrr::map2_dfc(
    pah_value_cols, lod_value_cols,
    \(pc, lc) !is.na(d[[pc]]) & !is.na(d[[lc]]) & d[[pc]] > d[[lc]]
  )

  d_tmp <- d |>
    dplyr::mutate(any_detect_row = rowSums(as.matrix(detect_mat), na.rm = TRUE) > 0)

  # collapse repeated measures first
  d_collapsed <- d_tmp |>
    dplyr::summarize(
      dplyr::across(dplyr::all_of(pah_value_cols), gm_nonzero_or_zero),
      any_detect = any(any_detect_row, na.rm = TRUE),
      .by = c(ParticipantID, Condition, Timing_new, SampleLocation)
    )

  # totalPAH from collapsed measured values
  d_run <- d_collapsed |>
    dplyr::mutate(
      n_nonmissing_pah = rowSums(!is.na(dplyr::pick(dplyr::all_of(pah_value_cols)))),
      totalPAH = rowSums(dplyr::pick(dplyr::all_of(pah_value_cols)), na.rm = TRUE),
      totalPAH = dplyr::if_else(n_nonmissing_pah == 0L, NA_real_, totalPAH)
    )

  c_const <- NA_real_

  if (response_mode == "all_data_plus_c") {
    min_pos <- suppressWarnings(min(d_run$totalPAH[d_run$totalPAH > 0], na.rm = TRUE))
    if (is.finite(min_pos)) c_const <- 0.5 * min_pos

    d_run <- d_run |>
      dplyr::mutate(
        response_log = dplyr::if_else(!is.na(totalPAH + c_const), log(totalPAH + c_const), NA_real_)
      )
  }

  if (response_mode == "detects_only") {
    d_run <- d_run |>
      dplyr::filter(any_detect, !is.na(totalPAH), totalPAH > 0) |>
      dplyr::mutate(response_log = log(totalPAH))
  }

  d_run <- d_run |>
    dplyr::mutate(
      run_id = run_id,
      run_label = run_label,
      pah_set = pah_set,
      response_mode = response_mode,
      c_const = c_const
    )

  run_diag <- tibble::tibble(
    run_id = run_id,
    run_label = run_label,
    pah_set = pah_set,
    response_mode = response_mode,
    n_rows_raw = nrow(d),
    n_rows_collapsed = nrow(d_collapsed),
    n_rows_model = nrow(d_run),
    n_participants = dplyr::n_distinct(d_run$ParticipantID),
    n_any_detect = sum(d_run$any_detect, na.rm = TRUE),
    n_response_nonmissing = sum(!is.na(d_run$response_log)),
    c_const = c_const
  )

  list(data = d_run, diag = run_diag, set_map = set_map)
}

run_objects <- purrr::pmap(
  list(
    analysis_grid$run_id,
    analysis_grid$run_label,
    analysis_grid$pah_set,
    analysis_grid$response_mode,
    analysis_grid$pah_names
  ),
  \(run_id, run_label, pah_set, response_mode, pah_names) {
    build_run_dataset(run_id, run_label, pah_set, response_mode, pah_names, d, detect_summary)
  }
)

analysis_data_list <- purrr::map(run_objects, "data")
names(analysis_data_list) <- analysis_grid$run_label

run_diagnostics <- purrr::map_dfr(run_objects, "diag") |>
  dplyr::arrange(run_id)

run_set_maps <- purrr::map(run_objects, "set_map")
names(run_set_maps) <- analysis_grid$run_label

step1c_artifacts <- list(
  analysis_data_list = analysis_data_list,
  run_diagnostics = run_diagnostics,
  run_set_maps = run_set_maps
)

saveRDS(step1c_artifacts, file.path(output_dir_10, "step1c_artifacts.rds"))

step1c_artifacts

# --------------------------------------------
# --- Step 1c (append): Tiny QC            ---
# --------------------------------------------

expected_collapsed_n <- d |>
  dplyr::distinct(ParticipantID, Condition, Timing_new, SampleLocation) |>
  nrow()

step1c_qc_tiny <- run_diagnostics |>
  dplyr::mutate(
    collapse_occurred = n_rows_collapsed < n_rows_raw,
    collapse_matches_expected = n_rows_collapsed == expected_collapsed_n,
    model_rows_valid = n_rows_model <= n_rows_collapsed
  ) |>
  dplyr::select(
    run_id, run_label, response_mode,
    n_rows_raw, n_rows_collapsed, n_rows_model,
    collapse_occurred, collapse_matches_expected, model_rows_valid
  )

stopifnot(all(step1c_qc_tiny$collapse_matches_expected))
stopifnot(all(step1c_qc_tiny$model_rows_valid))

saveRDS(step1c_qc_tiny, file.path(output_dir_10, "step1c_qc_tiny.rds"))

step1c_qc_tiny



# ---------------------------------
# --- Step 2: Fit baseline LMMs  ---
# ---         (revised)          ---
# ---------------------------------

prepare_model_df <- function(df) {
  model_vars <- c("response_log", "Condition", "Timing_new", "SampleLocation", "ParticipantID")

  df |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(model_vars), ~ !is.na(.x))) |>
    dplyr::mutate(
      ParticipantID = as.factor(ParticipantID),
      Condition = as.factor(Condition),
      Timing_new = as.factor(Timing_new),
      SampleLocation = as.factor(SampleLocation)
    ) |>
    droplevels()
}

run_lme_attempt <- function(df, use_varIdent = TRUE) {
  warn_msgs <- character(0)

  fit <- withCallingHandlers(
    tryCatch(
      {
        if (use_varIdent) {
          nlme::lme(
            fixed = response_log ~ Condition + Timing_new + SampleLocation,
            random = ~ 1 | ParticipantID,
            weights = nlme::varIdent(form = ~ 1 | SampleLocation),
            data = df,
            method = "REML",
            na.action = na.omit,
            control = nlme::lmeControl(opt = "optim")
          )
        } else {
          nlme::lme(
            fixed = response_log ~ Condition + Timing_new + SampleLocation,
            random = ~ 1 | ParticipantID,
            data = df,
            method = "REML",
            na.action = na.omit,
            control = nlme::lmeControl(opt = "optim")
          )
        }
      },
      error = \(e) e
    ),
    warning = \(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(fit, "error")) {
    list(model = NULL, error = conditionMessage(fit), warnings = warn_msgs)
  } else {
    list(model = fit, error = NA_character_, warnings = warn_msgs)
  }
}

fit_lme_run <- function(df, run_label) {
  df_model <- prepare_model_df(df)

  # pre-fit viability checks
  n_obs <- nrow(df_model)
  n_pid <- dplyr::n_distinct(df_model$ParticipantID)
  n_cond <- dplyr::n_distinct(df_model$Condition)
  n_time <- dplyr::n_distinct(df_model$Timing_new)
  n_loc <- dplyr::n_distinct(df_model$SampleLocation)

  viable <- n_obs > 0 && n_pid >= 2 && n_cond >= 2 && n_time >= 2 && n_loc >= 2

  if (!viable) {
    qc <- tibble::tibble(
      run_label = run_label,
      fit_success = FALSE,
      fit_path = "not_fit",
      variance_structure = NA_character_,
      fallback_used = FALSE,
      n_obs_model = n_obs,
      n_participant_model = n_pid,
      n_condition_model = n_cond,
      n_timing_model = n_time,
      n_location_model = n_loc,
      primary_error = "Pre-fit viability check failed (insufficient support).",
      fallback_error = NA_character_,
      primary_warnings = "",
      fallback_warnings = "",
      comparable_full = FALSE,
      comparable_any = FALSE,
      comparison_tier = "tier3_not_fit"
    )
    return(list(model = NULL, qc = qc))
  }

  # attempt 1: heteroscedastic
  attempt_hetero <- run_lme_attempt(df_model, use_varIdent = TRUE)

  if (!is.null(attempt_hetero$model)) {
    qc <- tibble::tibble(
      run_label = run_label,
      fit_success = TRUE,
      fit_path = "primary",
      variance_structure = "heteroscedastic_varIdent",
      fallback_used = FALSE,
      n_obs_model = n_obs,
      n_participant_model = n_pid,
      n_condition_model = n_cond,
      n_timing_model = n_time,
      n_location_model = n_loc,
      primary_error = NA_character_,
      fallback_error = NA_character_,
      primary_warnings = paste(attempt_hetero$warnings, collapse = " | "),
      fallback_warnings = "",
      comparable_full = TRUE,
      comparable_any = TRUE,
      comparison_tier = "tier1_full"
    )
    return(list(model = attempt_hetero$model, qc = qc))
  }

  # attempt 2: homoscedastic fallback
  attempt_homo <- run_lme_attempt(df_model, use_varIdent = FALSE)

  if (!is.null(attempt_homo$model)) {
    qc <- tibble::tibble(
      run_label = run_label,
      fit_success = TRUE,
      fit_path = "fallback",
      variance_structure = "homoscedastic",
      fallback_used = TRUE,
      n_obs_model = n_obs,
      n_participant_model = n_pid,
      n_condition_model = n_cond,
      n_timing_model = n_time,
      n_location_model = n_loc,
      primary_error = attempt_hetero$error,
      fallback_error = NA_character_,
      primary_warnings = paste(attempt_hetero$warnings, collapse = " | "),
      fallback_warnings = paste(attempt_homo$warnings, collapse = " | "),
      comparable_full = FALSE,
      comparable_any = TRUE,
      comparison_tier = "tier2_fallback"
    )
    return(list(model = attempt_homo$model, qc = qc))
  }

  qc <- tibble::tibble(
    run_label = run_label,
    fit_success = FALSE,
    fit_path = "failed_both",
    variance_structure = NA_character_,
    fallback_used = TRUE,
    n_obs_model = n_obs,
    n_participant_model = n_pid,
    n_condition_model = n_cond,
    n_timing_model = n_time,
    n_location_model = n_loc,
    primary_error = attempt_hetero$error,
    fallback_error = attempt_homo$error,
    primary_warnings = paste(attempt_hetero$warnings, collapse = " | "),
    fallback_warnings = paste(attempt_homo$warnings, collapse = " | "),
    comparable_full = FALSE,
    comparable_any = FALSE,
    comparison_tier = "tier3_not_fit"
  )

  list(model = NULL, qc = qc)
}

extract_fitstats <- function(mod, run_label) {
  tibble::tibble(
    run_label = run_label,
    n_obs = stats::nobs(mod),
    aic = AIC(mod),
    bic = BIC(mod),
    logLik = as.numeric(logLik(mod))
  )
}

extract_fixef <- function(mod, run_label) {
  tt <- summary(mod)$tTable
  tibble::tibble(
    run_label = run_label,
    term = rownames(tt),
    estimate = tt[, "Value"],
    std_error = tt[, "Std.Error"],
    df = tt[, "DF"],
    t_value = tt[, "t-value"],
    p_value = tt[, "p-value"]
  )
}

extract_condition_contrasts_safe <- function(mod, run_label) {
  tryCatch(
    {
      emmeans::emmeans(mod, "Condition") |>
        emmeans::contrast(method = "pairwise", adjust = "tukey") |>
        as.data.frame() |>
        tibble::as_tibble() |>
        dplyr::mutate(run_label = run_label, extract_error = NA_character_, .before = 1)
    },
    error = \(e) {
      tibble::tibble(
        run_label = run_label,
        contrast = NA_character_,
        estimate = NA_real_,
        SE = NA_real_,
        df = NA_real_,
        t.ratio = NA_real_,
        p.value = NA_real_,
        extract_error = conditionMessage(e)
      )
    }
  )
}

extract_timing_contrasts_safe <- function(mod, run_label) {
  tryCatch(
    {
      emmeans::emmeans(mod, "Timing_new") |>
        emmeans::contrast(method = "pairwise", adjust = "none") |>
        as.data.frame() |>
        tibble::as_tibble() |>
        dplyr::mutate(run_label = run_label, extract_error = NA_character_, .before = 1)
    },
    error = \(e) {
      tibble::tibble(
        run_label = run_label,
        contrast = NA_character_,
        estimate = NA_real_,
        SE = NA_real_,
        df = NA_real_,
        t.ratio = NA_real_,
        p.value = NA_real_,
        extract_error = conditionMessage(e)
      )
    }
  )
}

# ---- Fit all runs ----
fit_results_step2 <- purrr::imap(
  analysis_data_list,
  \(df, run_label) fit_lme_run(df, run_label)
)

models_step2 <- purrr::map(fit_results_step2, "model")
model_qc_step2 <- purrr::map_dfr(fit_results_step2, "qc")

models_step2_ok <- models_step2 |>
  purrr::keep(~ inherits(.x, "lme"))

# ---- Extract summaries from successful fits ----
fitstats_step2 <- purrr::imap_dfr(models_step2_ok, extract_fitstats)
fixef_step2 <- purrr::imap_dfr(models_step2_ok, extract_fixef)

condition_contrast_results <- purrr::imap_dfr(models_step2_ok, extract_condition_contrasts_safe) |>
  dplyr::mutate(significant = !is.na(p.value) & p.value < 0.05)

timing_results <- purrr::imap_dfr(models_step2_ok, extract_timing_contrasts_safe) |>
  dplyr::mutate(significant = !is.na(p.value) & p.value < 0.05)

# ---- Add run metadata + comparability lens ----
run_meta <- analysis_grid |>
  dplyr::select(run_id, run_label, pah_set, n_pahs, response_mode)

model_qc_step2 <- model_qc_step2 |>
  dplyr::left_join(run_meta, by = "run_label") |>
  dplyr::arrange(run_id)

fitstats_step2 <- fitstats_step2 |>
  dplyr::left_join(run_meta, by = "run_label") |>
  dplyr::left_join(model_qc_step2 |>
    dplyr::select(run_label, variance_structure, fit_path, comparison_tier, comparable_full, comparable_any),
    by = "run_label"
  ) |>
  dplyr::arrange(run_id)

fixef_step2 <- fixef_step2 |>
  dplyr::left_join(run_meta, by = "run_label") |>
  dplyr::left_join(model_qc_step2 |>
    dplyr::select(run_label, variance_structure, fit_path, comparison_tier, comparable_full, comparable_any),
    by = "run_label"
  ) |>
  dplyr::arrange(run_id, term)

condition_contrast_results <- condition_contrast_results |>
  dplyr::left_join(run_meta, by = "run_label") |>
  dplyr::left_join(model_qc_step2 |>
    dplyr::select(run_label, variance_structure, fit_path, comparison_tier, comparable_full, comparable_any),
    by = "run_label"
  ) |>
  dplyr::arrange(contrast, run_id)

timing_results <- timing_results |>
  dplyr::left_join(run_meta, by = "run_label") |>
  dplyr::left_join(model_qc_step2 |>
    dplyr::select(run_label, variance_structure, fit_path, comparison_tier, comparable_full, comparable_any),
    by = "run_label"
  ) |>
  dplyr::arrange(contrast, run_id)

# ---- Save Step 2 artifacts ----
step2_artifacts <- list(
  fit_results = fit_results_step2,
  models = models_step2,
  models_ok = models_step2_ok,
  model_qc = model_qc_step2,
  fitstats = fitstats_step2,
  fixef = fixef_step2,
  condition_contrast_results = condition_contrast_results,
  timing_results = timing_results
)

saveRDS(step2_artifacts, file.path(output_dir_10, "step2_artifacts.rds"))
saveRDS(model_qc_step2, file.path(output_dir_10, "step2_model_qc.rds"))

# Compact QC summary for immediate review
step2_qc <- model_qc_step2 |>
  dplyr::count(comparison_tier, variance_structure, fit_success, name = "n_runs")

saveRDS(step2_qc, file.path(output_dir_10, "step2_qc.rds"))

step2_qc

model_qc_step2 #yes, the all_data models required fallback to the homoscedastic structure
timing_results # doffing > donning in basically all cases
print(condition_contrast_results, n=Inf) # SS > OL and SS > SL across the board; unclear for SL vs. OL



# ---------------------------------
# --- Step 3: Boundary summaries ---
# ---------------------------------

# 3a) Prep contrast tables with CI and run order
run_order <- analysis_grid |>
  dplyr::arrange(response_mode, n_pahs) |>
  dplyr::pull(run_label)

condition_step3 <- condition_contrast_results |>
  dplyr::filter(!is.na(contrast), !is.na(estimate), !is.na(SE), !is.na(df)) |>
  dplyr::mutate(
    run_label = factor(run_label, levels = run_order),
    lower = estimate - stats::qt(0.975, df = df) * SE,
    upper = estimate + stats::qt(0.975, df = df) * SE
  )

timing_step3 <- timing_results |>
  dplyr::filter(!is.na(contrast), !is.na(estimate), !is.na(SE), !is.na(df)) |>
  dplyr::mutate(
    run_label = factor(run_label, levels = run_order),
    lower = estimate - stats::qt(0.975, df = df) * SE,
    upper = estimate + stats::qt(0.975, df = df) * SE
  )

# 3b) Boundary summary helper
summarize_boundary <- function(df) {
  df |>
    dplyr::group_by(comparison_tier, contrast) |>
    dplyr::summarize(
      n_runs = dplyr::n(),
      n_sig = sum(significant, na.rm = TRUE),
      prop_sig = n_sig / n_runs,
      est_min = min(estimate, na.rm = TRUE),
      est_max = max(estimate, na.rm = TRUE),
      abs_est_min = min(abs(estimate), na.rm = TRUE),
      abs_est_max = max(abs(estimate), na.rm = TRUE),
      all_positive = all(estimate > 0, na.rm = TRUE),
      all_negative = all(estimate < 0, na.rm = TRUE),
      direction_stable = all_positive | all_negative,
      .groups = "drop"
    ) |>
    dplyr::arrange(contrast, comparison_tier)
}

condition_boundary <- summarize_boundary(condition_step3)
timing_boundary <- summarize_boundary(timing_step3)

# 3c) Optional "comparable-only" boundary tables
condition_boundary_full <- condition_step3 |>
  dplyr::filter(comparable_full) |>
  summarize_boundary()

timing_boundary_full <- timing_step3 |>
  dplyr::filter(comparable_full) |>
  summarize_boundary()

# 3d) Visualization
p_condition_step3 <- condition_step3 |>
  ggplot2::ggplot(ggplot2::aes(x = run_label, y = estimate, color = significant)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.5) +
  ggplot2::facet_grid(contrast ~ comparison_tier, scales = "free_y") +
  ggplot2::coord_flip() +
  ggplot2::scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "#999999")) +
  ggplot2::labs(
    title = "Condition contrasts across boundary runs",
    x = NULL, y = "Estimate (log scale)", color = "p < 0.05"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

p_timing_step3 <- timing_step3 |>
  ggplot2::ggplot(ggplot2::aes(x = run_label, y = estimate, color = significant)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.5) +
  ggplot2::facet_grid(contrast ~ comparison_tier, scales = "free_y") +
  ggplot2::coord_flip() +
  ggplot2::scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#999999")) +
  ggplot2::labs(
    title = "Timing contrasts across boundary runs",
    x = NULL, y = "Estimate (log scale)", color = "p < 0.05"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

step3_plots <- p_condition_step3 / p_timing_step3

# 3e) Bundle artifacts
step3_artifacts <- list(
  condition_step3 = condition_step3,
  timing_step3 = timing_step3,
  condition_boundary = condition_boundary,
  timing_boundary = timing_boundary,
  condition_boundary_full = condition_boundary_full,
  timing_boundary_full = timing_boundary_full,
  p_condition_step3 = p_condition_step3,
  p_timing_step3 = p_timing_step3,
  step3_plots = step3_plots
)

saveRDS(step3_artifacts, file.path(output_dir_10, "step3_artifacts.rds"))

step3_artifacts


# ---------------------------------
# --- Step 4: Save outputs       ---
# ---------------------------------

# Core artifacts (.rds)
saveRDS(step0_setup, file.path(output_dir_10, "step0_setup.rds"))
saveRDS(step1_artifacts, file.path(output_dir_10, "step1_artifacts.rds"))
saveRDS(step1_qc, file.path(output_dir_10, "step1_qc.rds"))
saveRDS(step1c_artifacts, file.path(output_dir_10, "step1c_artifacts.rds"))
saveRDS(step1c_qc_tiny, file.path(output_dir_10, "step1c_qc_tiny.rds"))

saveRDS(step2_artifacts, file.path(output_dir_10, "step2_artifacts.rds"))
saveRDS(model_qc_step2, file.path(output_dir_10, "step2_model_qc.rds"))
saveRDS(step2_qc, file.path(output_dir_10, "step2_qc.rds"))

saveRDS(step3_artifacts, file.path(output_dir_10, "step3_artifacts.rds"))
saveRDS(condition_step3, file.path(output_dir_10, "condition_step3.rds"))
saveRDS(timing_step3, file.path(output_dir_10, "timing_step3.rds"))
saveRDS(condition_boundary, file.path(output_dir_10, "condition_boundary.rds"))
saveRDS(timing_boundary, file.path(output_dir_10, "timing_boundary.rds"))

# Plots (.png)
ggplot2::ggsave(
  filename = file.path(output_dir_10, "step3_condition_contrasts.png"),
  plot = p_condition_step3,
  width = 12, height = 8, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_10, "step3_timing_contrasts.png"),
  plot = p_timing_step3,
  width = 12, height = 5, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_10, "step3_combined.png"),
  plot = step3_plots,
  width = 14, height = 12, dpi = 300
)

# Optional manifest
step4_manifest <- tibble::tibble(
  file = list.files(output_dir_10, full.names = FALSE)
) |>
  dplyr::arrange(file)

saveRDS(step4_manifest, file.path(output_dir_10, "step4_manifest.rds"))

step4_manifest