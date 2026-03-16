# ============================================================
# 11_pah_sensitivity.R
# Explore sensitivity of results to the number of PAHs used in the total PAH metric
# Author: Nathan M. Drew
# Date: 2026-03-13
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


# ---------------------------------
# --- Step 1: Define PAH sets    ---
# ---------------------------------

# PAH/LOD columns from cleaned_data structure
pah_cols <- 8:22
lod_cols <- 24:38

stopifnot(length(pah_cols) == 15, length(lod_cols) == 15)
stopifnot(all(pah_cols %in% seq_along(d)), all(lod_cols %in% seq_along(d)))

pah_lookup_15 <- tibble(
  pah_name = names(d)[pah_cols],
  pah_col = names(d)[pah_cols],
  lod_col = names(d)[lod_cols]
)

# Global detect-rate summary (fixed once, based on raw measured data)
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
      pct_detect = dplyr::if_else(n_nonmissing > 0, 100 * n_detect / n_nonmissing, NA_real_)
    )
  }
) |>
  dplyr::arrange(dplyr::desc(pct_detect)) |>
  dplyr::mutate(rank_detect = dplyr::row_number())

# Freeze PAH sets: top1..top10 plus top15
set_sizes <- c(1:10, 15)

pah_sets <- purrr::set_names(
  set_sizes,
  nm = paste0("top", set_sizes)
) |>
  purrr::map(\(k) detect_summary |> dplyr::slice_head(n = k) |> dplyr::pull(pah_name))

pah_sets_map <- purrr::imap(
  pah_sets,
  \(set_names, set_id) {
    detect_summary |>
      dplyr::filter(pah_name %in% set_names) |>
      dplyr::mutate(set_id = set_id) |>
      dplyr::select(set_id, rank_detect, pah_name, pah_col, lod_col, pct_detect) |>
      dplyr::arrange(rank_detect)
  }
)

# Tiny QC
step1_qc <- tibble(
  set_id = names(pah_sets),
  expected_n = as.integer(stringr::str_remove(names(pah_sets), "top")),
  observed_n = purrr::map_int(pah_sets, length)
) |>
  dplyr::mutate(
    n_ok = expected_n == observed_n,
    top1_is_rank1 = pah_sets$top1[[1]] == detect_summary$pah_name[1]
  )

stopifnot(all(step1_qc$n_ok), all(step1_qc$top1_is_rank1))

step1_artifacts <- list(
  pah_lookup_15 = pah_lookup_15,
  detect_summary = detect_summary,
  pah_sets = pah_sets,
  pah_sets_map = pah_sets_map,
  step1_qc = step1_qc
)

step1_artifacts


# ---------------------------------
# --- Step 2: Build run grid     ---
# ---------------------------------

# Build one analysis run per PAH set
analysis_grid <- tibble(
  set_id = names(pah_sets),
  pah_names = unname(pah_sets)
) |>
  mutate(
    n_pahs = purrr::map_int(pah_names, length),
    set_num = as.integer(stringr::str_remove(set_id, "top")),
    run_id = row_number(),
    run_label = paste0("run_", run_id, "_", set_id),
    is_priority_set = set_id %in% c("top1", "top4", "top6", "top15")
  ) |>
  # attach PAH and LOD column vectors for each run
  mutate(
    pah_cols_run = purrr::map(pah_names, \(nm) {
      detect_summary |>
        filter(pah_name %in% nm) |>
        arrange(rank_detect) |>
        pull(pah_col)
    }),
    lod_cols_run = purrr::map(pah_names, \(nm) {
      detect_summary |>
        filter(pah_name %in% nm) |>
        arrange(rank_detect) |>
        pull(lod_col)
    })
  ) |>
  arrange(set_num) |>
  select(
    run_id, run_label, set_id, set_num, n_pahs, is_priority_set,
    pah_names, pah_cols_run, lod_cols_run
  )

# Tiny QC
step2_qc <- analysis_grid |>
  transmute(
    run_label,
    set_id,
    n_pahs_expected = set_num,
    n_pahs_observed = n_pahs,
    n_pah_cols = purrr::map_int(pah_cols_run, length),
    n_lod_cols = purrr::map_int(lod_cols_run, length),
    counts_ok = (n_pahs_observed == n_pahs_expected) &
      (n_pah_cols == n_pahs_observed) &
      (n_lod_cols == n_pahs_observed)
  )

stopifnot(
  nrow(analysis_grid) == 11,
  all(step2_qc$counts_ok),
  all(sort(analysis_grid$set_num) == c(1:10, 15))
)

step2_artifacts <- list(
  analysis_grid = analysis_grid,
  step2_qc = step2_qc
)

step2_artifacts



# ---------------------------------
# --- Step 3: MLE imputation     ---
# ---         (EnvStats)         ---
# ---------------------------------

# Helper: PAH-level beta from censored lognormal fit at reference LOD
calc_beta_from_fit <- function(meanlog, sdlog, lod_ref) {
  z1 <- (log(lod_ref) - meanlog) / sdlog
  z2 <- (log(lod_ref) - meanlog - sdlog^2) / sdlog
  ex_trunc <- exp(meanlog + 0.5 * sdlog^2) * stats::pnorm(z2) / stats::pnorm(z1)
  ex_trunc / lod_ref
}

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
    # no detects: fallback
    x_imp[ok] <- 0.5 * lod_ok
    return(list(
      x_imp = x_imp, beta = 0.5, meanlog = meanlog, sdlog = sdlog,
      fit_method = "all_censored_fallback_half_lod", fit_ok = FALSE,
      n_total = n_total, n_detect = n_detect, n_cens = n_cens
    ))
  }

  if (n_cens == 0) {
    # no censored values: no imputation needed
    return(list(
      x_imp = x_imp, beta = NA_real_, meanlog = NA_real_, sdlog = NA_real_,
      fit_method = "no_censoring", fit_ok = TRUE,
      n_total = n_total, n_detect = n_detect, n_cens = n_cens
    ))
  }

  # Censored lognormal fit (EnvStats)
  # Use censoring limits for censored observations (script 04 style)
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
    x_imp[ok][censored] <- 0.5 * lod_ok[censored]
    return(list(
      x_imp = x_imp, beta = 0.5, meanlog = meanlog, sdlog = sdlog,
      fit_method = "fit_failed_fallback_half_lod", fit_ok = FALSE,
      n_total = n_total, n_detect = n_detect, n_cens = n_cens
    ))
  }

  # Extract fitted parameters
  meanlog <- unname(fit$parameters[["meanlog"]])
  sdlog <- unname(fit$parameters[["sdlog"]])

  # PAH-level beta using reference LOD (median of all non-missing LODs; script 04 style)
  lod_ref <- stats::median(lod_ok, na.rm = TRUE)
  if (!is.finite(lod_ref) || lod_ref <= 0) {
    lod_ref <- stats::median(lod_ok[lod_ok > 0], na.rm = TRUE)
  }
  beta <- calc_beta_from_fit(meanlog = meanlog, sdlog = sdlog, lod_ref = lod_ref)

  # Impute only censored rows
  x_ok_imp <- x_ok
  x_ok_imp[censored] <- beta * lod_ok[censored]
  x_imp[ok] <- x_ok_imp

  list(
    x_imp = x_imp, beta = beta, meanlog = meanlog, sdlog = sdlog,
    fit_method = "elnormCensored_mle_beta_times_lod", fit_ok = TRUE,
    n_total = n_total, n_detect = n_detect, n_cens = n_cens
  )
}

# Build one imputed dataset per run (top1..top10, top15)
imputed_run_objects <- purrr::pmap(
  list(
    analysis_grid$run_label,
    analysis_grid$set_id,
    analysis_grid$pah_cols_run,
    analysis_grid$lod_cols_run
  ),
  \(run_label, set_id, pah_cols_run, lod_cols_run) {
    d_run <- d
    imp_log <- vector("list", length(pah_cols_run))
    imp_cols <- paste0(pah_cols_run, "_imp")

    for (i in seq_along(pah_cols_run)) {
      pc <- pah_cols_run[i]
      lc <- lod_cols_run[i]
      ic <- imp_cols[i]

      res <- impute_one_pah_envstats(d_run[[pc]], d_run[[lc]])
      d_run[[ic]] <- res$x_imp

      imp_log[[i]] <- tibble(
        run_label = run_label,
        set_id = set_id,
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

    d_run <- d_run |>
      dplyr::mutate(
        totalPAH_imputed = rowSums(dplyr::pick(dplyr::all_of(imp_cols)), na.rm = TRUE),
        totalPAH_imputed = dplyr::if_else(
          rowSums(!is.na(dplyr::pick(dplyr::all_of(imp_cols)))) == 0L,
          NA_real_,
          totalPAH_imputed
        ),
        run_label = run_label,
        set_id = set_id
      )

    list(
      data = d_run,
      imputation_log = dplyr::bind_rows(imp_log),
      imp_cols = imp_cols
    )
  }
)

imputed_data_list <- purrr::map(imputed_run_objects, "data")
names(imputed_data_list) <- analysis_grid$run_label

imputation_log <- purrr::map_dfr(imputed_run_objects, "imputation_log")

# Tiny QC
step3_qc <- imputation_log |>
  dplyr::group_by(run_label, set_id) |>
  dplyr::summarize(
    n_pahs = dplyr::n(),
    n_fit_ok = sum(fit_ok, na.rm = TRUE),
    any_fallback = any(!fit_ok),
    .groups = "drop"
  ) |>
  dplyr::arrange(run_label)

step3_artifacts <- list(
  imputed_data_list = imputed_data_list,
  imputation_log = imputation_log,
  step3_qc = step3_qc
)

step3_artifacts



# --------------------------------------------
# --- Step 2a (append): Top-6 beta QC       ---
# --------------------------------------------

top6_beta_prev <- readRDS("04_output/top6_beta.rds")

# Identify previous beta column safely
prev_beta_col <- if ("beta_mle" %in% names(top6_beta_prev)) {
  "beta_mle"
} else if ("beta" %in% names(top6_beta_prev)) {
  "beta"
} else {
  stop("No beta column found in 04_output/top6_beta.rds (expected 'beta_mle' or 'beta').")
}

# Recompute current top-6 betas using same EnvStats logic
calc_beta_envstats <- function(x, lod) {
  ok <- !is.na(x) & !is.na(lod)
  x_ok <- x[ok]
  lod_ok <- lod[ok]

  censored <- x_ok <= lod_ok
  n_detect <- sum(!censored)
  n_cens <- sum(censored)

  if (length(x_ok) == 0 || n_detect == 0 || n_cens == 0) return(NA_real_)

  # EnvStats expects censoring limits for censored observations
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
  if (inherits(fit, "try-error")) return(NA_real_)

  meanlog <- as.numeric(fit$parameters["meanlog"])
  sdlog <- as.numeric(fit$parameters["sdlog"])
  if (!is.finite(meanlog) || !is.finite(sdlog) || sdlog <= 0) return(NA_real_)

  lod_ref <- stats::median(lod_ok, na.rm = TRUE)
  if (!is.finite(lod_ref) || lod_ref <= 0) {
    lod_ref <- stats::median(lod_ok[lod_ok > 0], na.rm = TRUE)
  }
  if (!is.finite(lod_ref) || lod_ref <= 0) return(NA_real_)

  z1 <- (log(lod_ref) - meanlog) / sdlog
  z2 <- (log(lod_ref) - meanlog - sdlog^2) / sdlog
  ex_trunc <- exp(meanlog + 0.5 * sdlog^2) * stats::pnorm(z2) / stats::pnorm(z1)

  ex_trunc / lod_ref
}

top6_now <- detect_summary |>
  dplyr::arrange(rank_detect) |>
  dplyr::slice_head(n = 6) |>
  dplyr::mutate(
    beta_mle_current = purrr::map2_dbl(
      pah_col, lod_col,
      \(pc, lc) calc_beta_envstats(d[[pc]], d[[lc]])
    ),
    pah_name_key = stringr::str_to_lower(stringr::str_squish(as.character(pah_name)))
  ) |>
  dplyr::select(pah_name, pah_name_key, pah_col, lod_col, beta_mle_current)

# Standardize previous table (join by normalized PAH name key)
top6_prev_std <- top6_beta_prev |>
  dplyr::mutate(
    pah_name = as.character(pah_name),
    pah_name_key = stringr::str_to_lower(stringr::str_squish(pah_name)),
    beta_mle_prev = as.numeric(.data[[prev_beta_col]])
  ) |>
  dplyr::select(pah_name, pah_name_key, beta_mle_prev)

qc_top6_beta_compare <- top6_now |>
  dplyr::left_join(
    top6_prev_std |> dplyr::select(pah_name_key, beta_mle_prev),
    by = "pah_name_key"
  ) |>
  dplyr::mutate(
    abs_diff = abs(beta_mle_current - beta_mle_prev),
    pct_diff = 100 * (beta_mle_current - beta_mle_prev) / beta_mle_prev
  ) |>
  dplyr::select(
    pah_name, pah_col, lod_col,
    beta_mle_current, beta_mle_prev, abs_diff, pct_diff
  )

# Safe summary (no crash if no complete pairs)
n_complete <- sum(stats::complete.cases(qc_top6_beta_compare$beta_mle_current, qc_top6_beta_compare$beta_mle_prev))

qc_top6_beta_summary <- tibble::tibble(
  n_top6 = nrow(qc_top6_beta_compare),
  n_compared = n_complete,
  cor_beta = if (n_complete >= 2) {
    stats::cor(qc_top6_beta_compare$beta_mle_current, qc_top6_beta_compare$beta_mle_prev, use = "complete.obs")
  } else {
    NA_real_
  },
  max_abs_diff = if (n_complete >= 1) max(qc_top6_beta_compare$abs_diff, na.rm = TRUE) else NA_real_,
  median_abs_diff = if (n_complete >= 1) stats::median(qc_top6_beta_compare$abs_diff, na.rm = TRUE) else NA_real_
)

# Optional unmatched check
qc_top6_unmatched <- qc_top6_beta_compare |>
  dplyr::filter(is.na(beta_mle_prev)) |>
  dplyr::select(pah_name, pah_col, lod_col)

step2a_qc <- list(
  compare = qc_top6_beta_compare,
  summary = qc_top6_beta_summary,
  unmatched = qc_top6_unmatched
)

step2a_qc


# --------------------------------------------
# --- Step 2a diagnostic: why current NA?   ---
# --------------------------------------------

top6_diag <- detect_summary |>
  dplyr::arrange(rank_detect) |>
  dplyr::slice_head(n = 6) |>
  dplyr::mutate(
    n_complete = purrr::map2_int(
      pah_col, lod_col,
      \(pc, lc) sum(!is.na(d[[pc]]) & !is.na(d[[lc]]))
    ),
    n_detect = purrr::map2_int(
      pah_col, lod_col,
      \(pc, lc) sum(!is.na(d[[pc]]) & !is.na(d[[lc]]) & d[[pc]] > d[[lc]])
    ),
    n_cens = purrr::map2_int(
      pah_col, lod_col,
      \(pc, lc) sum(!is.na(d[[pc]]) & !is.na(d[[lc]]) & d[[pc]] <= d[[lc]])
    )
  ) |>
  dplyr::mutate(
    can_fit_mle = n_detect > 0 & n_cens > 0
  )

top6_diag



# ---------------------------------
# --- Step 4: GM collapse        ---
# ---------------------------------

gm_safe <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) NA_real_ else exp(mean(log(x)))
}

collapse_one_run <- function(df, run_label, set_id, pah_cols_run) {
  imp_cols <- paste0(pah_cols_run, "_imp")

  d_collapsed <- df |>
    dplyr::summarize(
      dplyr::across(dplyr::all_of(imp_cols), gm_safe),
      .by = c(ParticipantID, Condition, Timing_new, SampleLocation)
    ) |>
    dplyr::mutate(
      n_nonmissing_imp = rowSums(!is.na(dplyr::pick(dplyr::all_of(imp_cols)))),
      totalPAH_imputed = rowSums(dplyr::pick(dplyr::all_of(imp_cols)), na.rm = TRUE),
      totalPAH_imputed = dplyr::if_else(n_nonmissing_imp == 0L, NA_real_, totalPAH_imputed),
      log_totalPAH_imp = dplyr::if_else(totalPAH_imputed > 0, log(totalPAH_imputed), NA_real_),
      run_label = run_label,
      set_id = set_id
    )

  diag <- tibble::tibble(
    run_label = run_label,
    set_id = set_id,
    n_rows_raw = nrow(df),
    n_rows_collapsed = nrow(d_collapsed),
    n_response_nonmissing = sum(!is.na(d_collapsed$log_totalPAH_imp))
  )

  list(data = d_collapsed, diag = diag)
}

collapsed_run_objects <- purrr::pmap(
  list(
    analysis_grid$run_label,
    analysis_grid$set_id,
    analysis_grid$pah_cols_run
  ),
  \(run_label, set_id, pah_cols_run) {
    collapse_one_run(
      df = imputed_data_list[[run_label]],
      run_label = run_label,
      set_id = set_id,
      pah_cols_run = pah_cols_run
    )
  }
)

imputed_data_collapsed <- purrr::map(collapsed_run_objects, "data")
names(imputed_data_collapsed) <- analysis_grid$run_label

step4_qc <- purrr::map_dfr(collapsed_run_objects, "diag") |>
  dplyr::arrange(run_label)

step4_artifacts <- list(
  imputed_data_collapsed = imputed_data_collapsed,
  step4_qc = step4_qc
)

step4_artifacts



# ---------------------------------
# --- Step 5: Fit mixed models   ---
# ---------------------------------

prepare_model_df <- function(df) {
  df |>
    dplyr::filter(!is.na(log_totalPAH_imp)) |>
    dplyr::mutate(
      ParticipantID = as.factor(ParticipantID),
      Condition = as.factor(Condition),
      Timing_new = as.factor(Timing_new),
      SampleLocation = as.factor(SampleLocation)
    ) |>
    droplevels()
}

fit_lme_attempt <- function(df, hetero = TRUE) {
  warn <- character(0)

  fit <- withCallingHandlers(
    tryCatch(
      {
        if (hetero) {
          nlme::lme(
            fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
            random = ~ 1 | ParticipantID,
            weights = nlme::varIdent(form = ~ 1 | SampleLocation),
            data = df,
            method = "REML",
            na.action = na.omit,
            control = nlme::lmeControl(opt = "optim")
          )
        } else {
          nlme::lme(
            fixed = log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
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
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(fit, "error")) {
    list(model = NULL, error = conditionMessage(fit), warnings = warn)
  } else {
    list(model = fit, error = NA_character_, warnings = warn)
  }
}

fit_one_run <- function(df, run_label) {
  dm <- prepare_model_df(df)

  a1 <- fit_lme_attempt(dm, hetero = TRUE)
  if (!is.null(a1$model)) {
    return(list(
      model = a1$model,
      qc = tibble::tibble(
        run_label = run_label,
        fit_success = TRUE,
        fit_path = "primary",
        variance_structure = "heteroscedastic_varIdent",
        primary_error = NA_character_,
        fallback_error = NA_character_,
        primary_warnings = paste(a1$warnings, collapse = " | "),
        fallback_warnings = "",
        n_obs_model = nrow(dm)
      )
    ))
  }

  a2 <- fit_lme_attempt(dm, hetero = FALSE)
  if (!is.null(a2$model)) {
    return(list(
      model = a2$model,
      qc = tibble::tibble(
        run_label = run_label,
        fit_success = TRUE,
        fit_path = "fallback",
        variance_structure = "homoscedastic",
        primary_error = a1$error,
        fallback_error = NA_character_,
        primary_warnings = paste(a1$warnings, collapse = " | "),
        fallback_warnings = paste(a2$warnings, collapse = " | "),
        n_obs_model = nrow(dm)
      )
    ))
  }

  list(
    model = NULL,
    qc = tibble::tibble(
      run_label = run_label,
      fit_success = FALSE,
      fit_path = "failed_both",
      variance_structure = NA_character_,
      primary_error = a1$error,
      fallback_error = a2$error,
      primary_warnings = paste(a1$warnings, collapse = " | "),
      fallback_warnings = paste(a2$warnings, collapse = " | "),
      n_obs_model = nrow(dm)
    )
  )
}

# Fit all runs
fit_objects <- purrr::imap(imputed_data_collapsed, fit_one_run)

models_step5 <- purrr::map(fit_objects, "model")
models_step5_ok <- models_step5 |> purrr::keep(~ inherits(.x, "lme"))
model_qc_step5 <- purrr::map_dfr(fit_objects, "qc") |>
  dplyr::left_join(
    analysis_grid |> dplyr::select(run_id, run_label, set_id, set_num, n_pahs, is_priority_set),
    by = "run_label"
  ) |>
  dplyr::arrange(run_id)

# Tiny QC summary
step5_qc <- model_qc_step5 |>
  dplyr::count(fit_success, variance_structure, fit_path, name = "n_runs")

step5_artifacts <- list(
  models = models_step5,
  models_ok = models_step5_ok,
  model_qc = model_qc_step5,
  step5_qc = step5_qc
)

step5_artifacts



# ---------------------------------
# --- Step 6: Extract summaries  ---
# ---------------------------------

extract_fitstats <- function(mod, run_label) {
  tibble(
    run_label = run_label,
    n_obs = stats::nobs(mod),
    aic = AIC(mod),
    bic = BIC(mod),
    logLik = as.numeric(logLik(mod))
  )
}

extract_fixef <- function(mod, run_label) {
  tt <- summary(mod)$tTable

  tibble(
    run_label = run_label,
    term = rownames(tt),
    estimate = tt[, "Value"],
    std_error = tt[, "Std.Error"],
    df = tt[, "DF"],
    t_value = tt[, "t-value"],
    p_value = tt[, "p-value"]
  ) |>
    mutate(
      fold_change = exp(estimate),
      pct_change = 100 * (exp(estimate) - 1)
    )
}

extract_condition_contrasts <- function(mod, run_label) {
  emmeans::emmeans(mod, "Condition") |>
    emmeans::contrast(method = "pairwise", adjust = "tukey") |>
    as.data.frame() |>
    tibble::as_tibble() |>
    mutate(
      run_label = run_label,
      significant = p.value < 0.05,
      .before = 1
    )
}

extract_timing_contrasts <- function(mod, run_label) {
  emmeans::emmeans(mod, "Timing_new") |>
    emmeans::contrast(method = "pairwise", adjust = "none") |>
    as.data.frame() |>
    tibble::as_tibble() |>
    mutate(
      run_label = run_label,
      significant = p.value < 0.05,
      .before = 1
    )
}

run_meta <- analysis_grid |>
  select(run_id, run_label, set_id, set_num, n_pahs, is_priority_set)

fitstats_step6 <- purrr::imap_dfr(models_step5_ok, extract_fitstats) |>
  left_join(run_meta, by = "run_label") |>
  left_join(model_qc_step5 |> select(run_label, variance_structure, fit_path), by = "run_label") |>
  arrange(run_id)

fixef_step6 <- purrr::imap_dfr(models_step5_ok, extract_fixef) |>
  left_join(run_meta, by = "run_label") |>
  left_join(model_qc_step5 |> select(run_label, variance_structure, fit_path), by = "run_label") |>
  arrange(run_id, term)

condition_contrast_results <- purrr::imap_dfr(models_step5_ok, extract_condition_contrasts) |>
  left_join(run_meta, by = "run_label") |>
  left_join(model_qc_step5 |> select(run_label, variance_structure, fit_path), by = "run_label") |>
  arrange(contrast, run_id)

timing_results <- purrr::imap_dfr(models_step5_ok, extract_timing_contrasts) |>
  left_join(run_meta, by = "run_label") |>
  left_join(model_qc_step5 |> select(run_label, variance_structure, fit_path), by = "run_label") |>
  arrange(contrast, run_id)

step6_artifacts <- list(
  fitstats = fitstats_step6,
  fixef = fixef_step6,
  condition_contrast_results = condition_contrast_results,
  timing_results = timing_results
)

# Tiny QC
step6_qc <- list(
  n_models = length(models_step5_ok),
  n_fitstats_rows = nrow(fitstats_step6),
  n_condition_rows = nrow(condition_contrast_results),
  n_timing_rows = nrow(timing_results)
)

step6_artifacts



# ---------------------------------
# --- Step 7: Boundary summaries ---
# ---------------------------------

# 7a) Condition boundary summary across all runs
condition_boundary <- condition_contrast_results |>
  dplyr::group_by(contrast) |>
  dplyr::summarize(
    n_runs = dplyr::n(),
    n_sig = sum(significant, na.rm = TRUE),
    prop_sig = n_sig / n_runs,
    est_min = min(estimate, na.rm = TRUE),
    est_max = max(estimate, na.rm = TRUE),
    all_positive = all(estimate > 0, na.rm = TRUE),
    all_negative = all(estimate < 0, na.rm = TRUE),
    direction_stable = all_positive | all_negative,
    .groups = "drop"
  ) |>
  dplyr::arrange(contrast)

# 7b) Timing boundary summary across all runs
timing_boundary <- timing_results |>
  dplyr::group_by(contrast) |>
  dplyr::summarize(
    n_runs = dplyr::n(),
    n_sig = sum(significant, na.rm = TRUE),
    prop_sig = n_sig / n_runs,
    est_min = min(estimate, na.rm = TRUE),
    est_max = max(estimate, na.rm = TRUE),
    all_positive = all(estimate > 0, na.rm = TRUE),
    all_negative = all(estimate < 0, na.rm = TRUE),
    direction_stable = all_positive | all_negative,
    .groups = "drop"
  )

# 7c) Priority-only summaries (top1, top4, top6, top15)
condition_boundary_priority <- condition_contrast_results |>
  dplyr::filter(is_priority_set) |>
  dplyr::group_by(contrast) |>
  dplyr::summarize(
    n_runs = dplyr::n(),
    n_sig = sum(significant, na.rm = TRUE),
    prop_sig = n_sig / n_runs,
    est_min = min(estimate, na.rm = TRUE),
    est_max = max(estimate, na.rm = TRUE),
    all_positive = all(estimate > 0, na.rm = TRUE),
    all_negative = all(estimate < 0, na.rm = TRUE),
    direction_stable = all_positive | all_negative,
    .groups = "drop"
  ) |>
  dplyr::arrange(contrast)

timing_boundary_priority <- timing_results |>
  dplyr::filter(is_priority_set) |>
  dplyr::group_by(contrast) |>
  dplyr::summarize(
    n_runs = dplyr::n(),
    n_sig = sum(significant, na.rm = TRUE),
    prop_sig = n_sig / n_runs,
    est_min = min(estimate, na.rm = TRUE),
    est_max = max(estimate, na.rm = TRUE),
    all_positive = all(estimate > 0, na.rm = TRUE),
    all_negative = all(estimate < 0, na.rm = TRUE),
    direction_stable = all_positive | all_negative,
    .groups = "drop"
  )

# 7d) Optional compact run-level table for quick review
condition_run_table <- condition_contrast_results |>
  dplyr::select(run_label, set_id, set_num, contrast, estimate, SE, p.value, significant, is_priority_set) |>
  dplyr::arrange(contrast, set_num)

timing_run_table <- timing_results |>
  dplyr::select(run_label, set_id, set_num, contrast, estimate, SE, p.value, significant, is_priority_set) |>
  dplyr::arrange(contrast, set_num)

step7_artifacts <- list(
  condition_boundary = condition_boundary,
  timing_boundary = timing_boundary,
  condition_boundary_priority = condition_boundary_priority,
  timing_boundary_priority = timing_boundary_priority,
  condition_run_table = condition_run_table,
  timing_run_table = timing_run_table
)

step7_artifacts



# ---------------------------------
# --- Step 8: Sensitivity plots  ---
# ---------------------------------

# Run order for plotting
set_order <- analysis_grid |>
  dplyr::arrange(set_num) |>
  dplyr::pull(set_id)

condition_plot_data <- condition_contrast_results |>
  dplyr::mutate(
    set_id = factor(set_id, levels = set_order),
    ci_low = estimate - 1.96 * SE,
    ci_high = estimate + 1.96 * SE
  )

timing_plot_data <- timing_results |>
  dplyr::mutate(
    set_id = factor(set_id, levels = set_order),
    ci_low = estimate - 1.96 * SE,
    ci_high = estimate + 1.96 * SE
  )

p_condition_step8 <- condition_plot_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = set_num, y = estimate,
      ymin = ci_low, ymax = ci_high,
      color = is_priority_set
    )
  ) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_pointrange(size = 0.4) +
  ggplot2::geom_line(ggplot2::aes(group = contrast), color = "grey70", linewidth = 0.4) +
  ggplot2::facet_wrap(~ contrast, scales = "free_y") +
  ggplot2::scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey60")) +
  ggplot2::scale_x_continuous(breaks = c(1:10, 15)) +
  ggplot2::labs(
    title = "Condition contrasts by PAH set size",
    x = "Top-k PAHs included in totalPAH_imputed",
    y = "Estimate (log scale)",
    color = "Priority set"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

p_timing_step8 <- timing_plot_data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = set_num, y = estimate,
      ymin = ci_low, ymax = ci_high,
      color = is_priority_set
    )
  ) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggplot2::geom_pointrange(size = 0.5) +
  ggplot2::geom_line(color = "grey70", linewidth = 0.5) +
  ggplot2::scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "grey60")) +
  ggplot2::scale_x_continuous(breaks = c(1:10, 15)) +
  ggplot2::labs(
    title = "Timing contrast by PAH set size",
    x = "Top-k PAHs included in totalPAH_imputed",
    y = "Estimate (log scale)",
    color = "Priority set"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

p_step8_combined <- p_condition_step8 / p_timing_step8 +
  patchwork::plot_annotation(title = "Sensitivity to PAH set size")

step8_artifacts <- list(
  p_condition_step8 = p_condition_step8,
  p_timing_step8 = p_timing_step8,
  p_step8_combined = p_step8_combined
)

step8_artifacts



# ---------------------------------
# --- Step 9: Save artifacts     ---
# ---------------------------------

output_dir_11 <- "11_output"
dir.create(output_dir_11, showWarnings = FALSE, recursive = TRUE)

# Save core artifacts (.rds)
saveRDS(step1_artifacts, file.path(output_dir_11, "step1_artifacts.rds"))
saveRDS(step2_artifacts, file.path(output_dir_11, "step2_artifacts.rds"))
saveRDS(step2a_qc, file.path(output_dir_11, "step2a_qc.rds"))
saveRDS(step3_artifacts, file.path(output_dir_11, "step3_artifacts.rds"))
saveRDS(step4_artifacts, file.path(output_dir_11, "step4_artifacts.rds"))
saveRDS(step5_artifacts, file.path(output_dir_11, "step5_artifacts.rds"))
saveRDS(step6_artifacts, file.path(output_dir_11, "step6_artifacts.rds"))
saveRDS(step7_artifacts, file.path(output_dir_11, "step7_artifacts.rds"))
saveRDS(step8_artifacts, file.path(output_dir_11, "step8_artifacts.rds"))

# Save commonly reused objects (.rds)
saveRDS(analysis_grid, file.path(output_dir_11, "analysis_grid.rds"))
saveRDS(detect_summary, file.path(output_dir_11, "detect_summary.rds"))
saveRDS(imputation_log, file.path(output_dir_11, "imputation_log.rds"))
saveRDS(step4_qc, file.path(output_dir_11, "step4_qc.rds"))
saveRDS(model_qc_step5, file.path(output_dir_11, "model_qc_step5.rds"))
saveRDS(condition_contrast_results, file.path(output_dir_11, "condition_contrast_results.rds"))
saveRDS(timing_results, file.path(output_dir_11, "timing_results.rds"))
saveRDS(condition_boundary, file.path(output_dir_11, "condition_boundary.rds"))
saveRDS(timing_boundary, file.path(output_dir_11, "timing_boundary.rds"))

# Save plots (.png)
ggplot2::ggsave(
  filename = file.path(output_dir_11, "step8_condition_contrasts.png"),
  plot = p_condition_step8,
  width = 11, height = 7, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_11, "step8_timing_contrast.png"),
  plot = p_timing_step8,
  width = 11, height = 5, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir_11, "step8_combined.png"),
  plot = p_step8_combined,
  width = 12, height = 10, dpi = 300
)

# Manifest
step9_manifest <- tibble::tibble(
  file = list.files(output_dir_11, full.names = FALSE)
) |>
  dplyr::arrange(file)

saveRDS(step9_manifest, file.path(output_dir_11, "step9_manifest.rds"))

step9_manifest