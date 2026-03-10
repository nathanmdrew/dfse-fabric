# ============================================================
# 05_beta_sensitivity_analysis.R
# Sensitivity analysis: effect of beta factor on model results
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

setwd("C:/Users/vom8/dfse-fabric/")

d_base <- readRDS("data/cleaned_data.RDS")

d_base <- d_base |>
  rename(SampleLocation = `Sample Location`) |>
  mutate(
    SampleLocation = factor(SampleLocation),
    Condition      = factor(Condition),
    Timing_new     = factor(Timing_new),
    ParticipantID  = factor(ParticipantID)
  )

# Collapse Shirt → Sleeve (PI-approved)
d_base <- d_base |>
  mutate(
    SampleLocation = fct_recode(SampleLocation, "Sleeve" = "Shirt") |> fct_drop()
  )

# Verify
d_base |> count(SampleLocation, sort = TRUE)
stopifnot(nrow(d_base) == 448)

# Rebuild pair_map and top 6
pah_cols <- 8:22
lod_cols <- 24:38

pair_map <- tibble(
  pair_id  = seq_along(pah_cols),
  pah_col  = pah_cols,
  lod_col  = lod_cols,
  pah_name = names(d_base)[pah_cols],
  lod_name = names(d_base)[lod_cols]
)

pah_detect_rank <- pair_map |>
  mutate(
    n_detect = map2_int(pah_col, lod_col, \(pc, lc) {
      sum(d_base[[pc]] > d_base[[lc]], na.rm = TRUE)
    })
  ) |>
  arrange(desc(n_detect))

top6 <- pah_detect_rank |> slice_head(n = 6)

# Load MLE betas from script 04
top6_beta <- readRDS("04_output/top6_beta.rds")

# Define beta grid
# Note: MLE betas are per-PAH (~0.11 to 0.21); conventional betas are single values
beta_grid <- tibble(
  beta_label = c("Near-zero (0.01)", "MLE (per-PAH)", "LOD/2 (0.5)", "LOD/√2 (0.707)", "Full LOD (1.0)"),
  beta_value = list(0.01, top6_beta$beta_mle, 0.5, 1/sqrt(2), 1.0),
  is_mle     = c(FALSE, TRUE, FALSE, FALSE, FALSE)
)

beta_grid

# -----------------------------------------------
# --- Step 1: Impute and fit for each beta    ---
# -----------------------------------------------

# Function: impute top 6 PAHs with a given beta, sum, log-transform, fit LME
fit_with_beta <- function(beta_vals, beta_label, is_mle, d, top6_beta) {

  # If single beta, replicate for all 6 PAHs
  if (length(beta_vals) == 1) {
    beta_vals <- rep(beta_vals, nrow(top6_beta))
  }

  # Impute
  for (i in seq_len(nrow(top6_beta))) {
    pc   <- top6_beta$pah_col[i]
    lc   <- top6_beta$lod_col[i]
    beta <- beta_vals[i]
    nm   <- top6_beta$pah_name[i]

    pah_vals <- d[[pc]]
    lod_vals <- d[[lc]]
    is_nd    <- !is.na(pah_vals) & !is.na(lod_vals) & (pah_vals <= lod_vals)

    imp_name <- paste0(nm, "_imp")
    d[[imp_name]] <- pah_vals
    d[[imp_name]][is_nd] <- beta * lod_vals[is_nd]
  }

  imp_cols <- paste0(top6_beta$pah_name, "_imp")
  d <- d |>
    mutate(
      totalPAH_imputed = rowSums(pick(all_of(imp_cols)), na.rm = TRUE),
      log_totalPAH_imp = log(totalPAH_imputed)
    )

  # Check for -Inf (from log(0)) — skip model if present
  if (any(is.infinite(d$log_totalPAH_imp))) {
    warning(paste0(beta_label, ": log(0) encountered, skipping model fit"))
    return(NULL)
  }

  # Fit LME with varIdent (same structure as m5)
  mod <- tryCatch(
    lme(
      log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
      random = ~1 | ParticipantID,
      weights = varIdent(form = ~1 | SampleLocation),
      data = d
    ),
    error = function(e) {
      warning(paste0(beta_label, ": model failed — ", e$message))
      NULL
    }
  )

  mod
}

# Fit models across the beta grid
models <- pmap(
  list(
    beta_vals  = beta_grid$beta_value,
    beta_label = beta_grid$beta_label,
    is_mle     = beta_grid$is_mle
  ),
  \(beta_vals, beta_label, is_mle)
    fit_with_beta(beta_vals, beta_label, is_mle, d_base, top6_beta)
)

names(models) <- beta_grid$beta_label

# Check which models converged
map_lgl(models, \(m) !is.null(m)) #Full LOD did not converge

# Refit Full LOD with homogeneous variance
models[["Full LOD (1.0)"]] <- fit_with_beta(
  beta_vals  = 1.0,
  beta_label = "Full LOD (1.0) [homog. var.]",
  is_mle     = FALSE,
  d          = d_base,
  top6_beta  = top6_beta
)

# But swap out the varIdent model for a simple lme
models[["Full LOD (1.0)"]] <- tryCatch(
  lme(
    log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
    random = ~1 | ParticipantID,
    data = {
      d_temp <- d_base
      for (i in seq_len(nrow(top6_beta))) {
        pc <- top6_beta$pah_col[i]; lc <- top6_beta$lod_col[i]
        nm <- top6_beta$pah_name[i]
        pah_vals <- d_temp[[pc]]; lod_vals <- d_temp[[lc]]
        is_nd <- !is.na(pah_vals) & !is.na(lod_vals) & (pah_vals <= lod_vals)
        d_temp[[paste0(nm, "_imp")]] <- pah_vals
        d_temp[[paste0(nm, "_imp")]][is_nd] <- 1.0 * lod_vals[is_nd]
      }
      imp_cols <- paste0(top6_beta$pah_name, "_imp")
      d_temp |>
        mutate(
          totalPAH_imputed = rowSums(pick(all_of(imp_cols)), na.rm = TRUE),
          log_totalPAH_imp = log(totalPAH_imputed)
        )
    }
  ),
  error = function(e) { warning(e$message); NULL }
)

is.null(models[["Full LOD (1.0)"]])

summary(models[["Full LOD (1.0)"]])
# So it actually did converge; however the variance structure is now homogeneous
# and not directly comparable to the other models.
# Results are also strange; for example, finger, fly, lower pant, etc. are negative WRT back
# May need to omit the full LOD.

# Update beta grid: replace Full LOD (1.0) with Near-LOD (0.99)
beta_grid <- tibble(
  beta_label = c("Near-zero (0.01)", "MLE (per-PAH)", "LOD/2 (0.5)", "LOD/√2 (0.707)", "Near-LOD (0.99)"),
  beta_value = list(0.01, top6_beta$beta_mle, 0.5, 1/sqrt(2), 0.99),
  is_mle     = c(FALSE, TRUE, FALSE, FALSE, FALSE)
)

# Refit just the Near-LOD model
models[["Near-LOD (0.99)"]] <- fit_with_beta(
  beta_vals  = 0.99,
  beta_label = "Near-LOD (0.99)",
  is_mle     = FALSE,
  d          = d_base,
  top6_beta  = top6_beta
)

# Remove the old Full LOD entry
models[["Full LOD (1.0)"]] <- NULL

# Check convergence
map_lgl(models, \(m) !is.null(m))

# Update Near-LOD to 0.9
beta_grid <- tibble(
  beta_label = c("Near-zero (0.01)", "MLE (per-PAH)", "LOD/2 (0.5)", "LOD/√2 (0.707)", "Near-LOD (0.9)"),
  beta_value = list(0.01, top6_beta$beta_mle, 0.5, 1/sqrt(2), 0.9),
  is_mle     = c(FALSE, TRUE, FALSE, FALSE, FALSE)
)

models[["Near-LOD (0.99)"]] <- NULL

models[["Near-LOD (0.9)"]] <- fit_with_beta(
  beta_vals  = 0.9,
  beta_label = "Near-LOD (0.9)",
  is_mle     = FALSE,
  d          = d_base,
  top6_beta  = top6_beta
)

map_lgl(models, \(m) !is.null(m))


# -----------------------------------------------
# --- Step 2: Compare fixed-effect estimates  ---
# ---         across beta values              ---
# -----------------------------------------------

# Extract Condition contrasts from each model
contrast_results <- imap_dfr(models, \(mod, label) {
  emm <- emmeans(mod, "Condition")
  contr <- contrast(emm, method = "pairwise", adjust = "tukey")
  contr_df <- as.data.frame(contr)
  contr_df$beta_label <- label
  contr_df
})

# Order beta labels from lowest to highest
beta_order <- c("Near-zero (0.01)", "MLE (per-PAH)", "LOD/2 (0.5)",
                "LOD/√2 (0.707)", "Near-LOD (0.9)")

contrast_results <- contrast_results |>
  mutate(
    beta_label = factor(beta_label, levels = beta_order),
    significant = p.value < 0.05
  )

# Summary table
contrast_results |>
  select(beta_label, contrast, estimate, SE, p.value, significant) |>
  arrange(contrast, beta_label)


# -----------------------------------------------
# --- Step 3: Visualize sensitivity results   ---
# -----------------------------------------------

# Condition contrasts: coefficient plot with CIs
p_condition <- contrast_results |>
  mutate(
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE
  ) |>
  ggplot(aes(x = beta_label, y = estimate, color = significant)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "#CC79A7")) +
  facet_wrap(~contrast, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Sensitivity of Condition Contrasts to Beta Factor",
    subtitle = "Estimates ± 95% CI on log scale (Tukey-adjusted p-values)",
    x = NULL, y = "Estimate (log scale)",
    color = "p < 0.05"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p_condition

# Also extract Timing contrasts across beta levels
timing_results <- imap_dfr(models, \(mod, label) {
  emm <- emmeans(mod, "Timing_new")
  contr <- contrast(emm, method = "pairwise", adjust = "none")
  contr_df <- as.data.frame(contr)
  contr_df$beta_label <- label
  contr_df
}) |>
  mutate(
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
    subtitle = "Donning/Firefighting vs. Doffing: Estimate ± 95% CI",
    x = NULL, y = "Estimate (log scale)",
    color = "p < 0.05"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p_condition / p_timing +
  plot_annotation(title = "Beta Sensitivity Analysis")


# -----------------------------------------------
# --- Step 4: Summarize and save artifacts    ---
# -----------------------------------------------

# Summary tables
saveRDS(contrast_results, "05_output/condition_contrasts_by_beta.rds")
saveRDS(timing_results, "05_output/timing_contrasts_by_beta.rds")
saveRDS(beta_grid, "05_output/beta_grid.rds")
saveRDS(models, "05_output/models_by_beta.rds")

# Plots
ggsave("05_output/sensitivity_condition_contrasts.png",
       p_condition, width = 12, height = 6, dpi = 300)

ggsave("05_output/sensitivity_timing_contrast.png",
       p_timing, width = 10, height = 5, dpi = 300)

ggsave("05_output/sensitivity_combined.png",
       p_condition / p_timing + plot_annotation(title = "Beta Sensitivity Analysis"),
       width = 12, height = 10, dpi = 300)

# Summary narrative table for the report
sensitivity_summary <- tribble(
  ~contrast,  ~direction,  ~robust,  ~notes,
  "SS vs OL", "SS > OL",   "Yes",   "Significant across all 5 beta values (p < 1e-05)",
  "SS vs SL", "SS > SL",   "Yes",   "Significant across all 5 beta values (p < 0.01)",
  "SL vs OL", "SL > OL",   "Partial", "Non-significant at near-zero and MLE betas; significant at LOD/2, LOD/√2, near-LOD",
  "Doffing vs Donning", "Doffing > Donning", "Yes", "Significant across all 5 beta values; magnitude sensitive to beta (larger effect at smaller beta)"
)

saveRDS(sensitivity_summary, "05_output/sensitivity_summary.rds")

sensitivity_summary