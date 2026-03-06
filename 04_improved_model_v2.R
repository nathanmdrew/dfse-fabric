### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-03-06
### Purpose: Improve on the imputed model by accounting for variance by location

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

# ---------------------------------
# --- Step 4: Fit the LMM       ---
# ---------------------------------

# No constant needed — all values are positive
d <- d |>
  mutate(log_totalPAH_imp = log(totalPAH_imputed))

m3 <- lme(log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
    random = ~1 | ParticipantID,
    weights = varIdent(form = ~1 | SampleLocation),
    data = d)

summary(m3)

# -----------------------------------------
# --- Step 5: Residual diagnostics (m3) ---
# -----------------------------------------

diag_df3 <- tibble(
  fitted = fitted(m3),
  resid  = residuals(m3, type = "normalized")
)

# 1. Residuals vs. fitted
p1 <- ggplot(diag_df3, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, color = "steelblue", linewidth = 0.8) +
  labs(title = "Residuals vs. Fitted", x = "Fitted values", y = "Normalized residuals") +
  theme_minimal(base_size = 12)

# 2. QQ plot of residuals
p2 <- ggplot(diag_df3, aes(sample = resid)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(color = "red") +
  labs(title = "Normal Q-Q", x = "Theoretical quantiles", y = "Normalized residuals") +
  theme_minimal(base_size = 12)

# 3. Histogram of residuals
p3 <- ggplot(diag_df3, aes(x = resid)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "black") +
  labs(title = "Residual Distribution", x = "Normalized residuals", y = "Count") +
  theme_minimal(base_size = 12)

# 4. Random effects QQ
re_df3 <- tibble(
  re = ranef(m3)[["(Intercept)"]]
)

p4 <- ggplot(re_df3, aes(sample = re)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(color = "red") +
  labs(title = "Random Effects Q-Q", x = "Theoretical quantiles",
       y = "Participant intercepts") +
  theme_minimal(base_size = 12)

# Combine
(p1 + p2) / (p3 + p4) +
  plot_annotation(title = "Model m3: Residual Diagnostics (Heterogeneous Variance)")

# ------------------------------------------------
# --- Step 6: Fit the Gamma GLMM with log link ---
# ------------------------------------------------

m4 <- glmmTMB(
  totalPAH_imputed ~ Condition + Timing_new + SampleLocation + (1 | ParticipantID),
  family = Gamma(link = "log"),
  dispformula = ~SampleLocation,
  data = d
)

summary(m4)

# Convergence warning
# Look at the Shirt dispersion estimate — 17.76 with SE of 37.3. 
# That's completely unidentifiable. You can't estimate a location-specific dispersion from 
# 1 observation. Lower Chest (n=2) is also suspect.

# -------------------------------------------------------
# --- Step 6b: Collapse ad hoc locations and refit m4 ---
# -------------------------------------------------------

# Collapse: Shirt → Sleeve, Lower Chest → Chest, Fly → Pant
d <- d |>
  mutate(
    SampleLocation = fct_recode(SampleLocation,
      "Sleeve" = "Shirt",
      "Chest"  = "Lower Chest",
      "Pant"   = "Fly"
    ) |> fct_drop()
  )

# Verify
d |> count(SampleLocation, sort = TRUE)

# Refit Gamma GLMM
m4b <- glmmTMB(
  totalPAH_imputed ~ Condition + Timing_new + SampleLocation + (1 | ParticipantID),
  family = Gamma(link = "log"),
  dispformula = ~SampleLocation,
  data = d
)

summary(m4b)

# ------------------------------------------------
# --- Step 7: Residual diagnostics (m4b)       ---
# ------------------------------------------------

# Pearson residuals are appropriate for Gamma GLMM
diag_df4b <- tibble(
  fitted = fitted(m4b),
  resid  = residuals(m4b, type = "pearson")
)

# 1. Residuals vs. fitted
p1 <- ggplot(diag_df4b, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, color = "steelblue", linewidth = 0.8) +
  labs(title = "Residuals vs. Fitted", x = "Fitted values", y = "Pearson residuals") +
  theme_minimal(base_size = 12)

# 2. QQ plot of deviance residuals (closer to normal for GLMMs)
diag_df4b_dev <- tibble(
  resid_dev = residuals(m4b, type = "deviance")
)

p2 <- ggplot(diag_df4b_dev, aes(sample = resid_dev)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(color = "red") +
  labs(title = "Normal Q-Q (Deviance Residuals)", x = "Theoretical quantiles",
       y = "Deviance residuals") +
  theme_minimal(base_size = 12)

# 3. Histogram of Pearson residuals
p3 <- ggplot(diag_df4b, aes(x = resid)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "black") +
  labs(title = "Residual Distribution", x = "Pearson residuals", y = "Count") +
  theme_minimal(base_size = 12)

# 4. Random effects QQ
re_df4b <- tibble(
  re = ranef(m4b)$cond$ParticipantID[["(Intercept)"]]
)

p4 <- ggplot(re_df4b, aes(sample = re)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(color = "red") +
  labs(title = "Random Effects Q-Q", x = "Theoretical quantiles",
       y = "Participant intercepts") +
  theme_minimal(base_size = 12)

# Combine
(p1 + p2) / (p3 + p4) +
  plot_annotation(title = "Model m4b: Residual Diagnostics (Gamma GLMM)")


# -----------------------------------------
# --- Step 8: Refit m3 with collapsed    ---
# ---         locations (m3b)            ---
# -----------------------------------------

m3b <- lme(
  log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
  random = ~1 | ParticipantID,
  weights = varIdent(form = ~1 | SampleLocation),
  data = d
)

summary(m3b)

# Some odd results: finger, lower pant, thumb are non-sig
# What do the Fly observations look like?
d_orig <- readRDS("data/cleaned_data.RDS")

d_orig |>
  filter(`Sample Location` == "Fly") |>
  select(`Sample Location`, Condition, Timing_new, ParticipantID, new_totalPAH)

# -----------------------------------------
# --- Step 9: Drop ad hoc locations,    ---
# ---         refit LME (m5)            ---
# -----------------------------------------

# Reload original data and rebuild from scratch
d <- readRDS("data/cleaned_data.RDS")

d <- d |>
  rename(SampleLocation = `Sample Location`) |>
  mutate(
    SampleLocation = factor(SampleLocation),
    Condition      = factor(Condition),
    Timing_new     = factor(Timing_new),
    ParticipantID  = factor(ParticipantID)
  )

# Collapse Shirt → Sleeve only
d <- d |>
  mutate(
    SampleLocation = fct_recode(SampleLocation, "Sleeve" = "Shirt") |> fct_drop()
  )

# Drop Lower Chest and Fly
d <- d |>
  filter(!SampleLocation %in% c("Lower Chest", "Fly")) |>
  mutate(SampleLocation = fct_drop(SampleLocation))

# Verify
d |> count(SampleLocation, sort = TRUE)

# Rebuild imputation (same as earlier steps)
for (i in seq_len(nrow(top6_beta))) {
  pc   <- top6_beta$pah_col[i]
  lc   <- top6_beta$lod_col[i]
  beta <- top6_beta$beta_mle[i]
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

# How many observations dropped?
cat("Observations:", nrow(d), "(dropped", 448 - nrow(d), ")\n")

# Fit m5: LME with heterogeneous variance, ad hoc locations removed
m5 <- lme(
  log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
  random = ~1 | ParticipantID,
  weights = varIdent(form = ~1 | SampleLocation),
  data = d
)

summary(m5)

# -----------------------------------------
# --- Step 10: Test for heterogeneous   ---
# ---          variance                 ---
# -----------------------------------------

# Fit homogeneous variance version with lme (same data as m5)
m5_homog <- lme(
  log_totalPAH_imp ~ Condition + Timing_new + SampleLocation,
  random = ~1 | ParticipantID,
  data = d
)

# Likelihood ratio test: homogeneous vs. heterogeneous
anova(m5_homog, m5)


# -----------------------------------------
# --- Step 11: Save artifacts to 04_output
# -----------------------------------------

# Models
saveRDS(m3, "04_output/m3_lme_varident.rds")
saveRDS(m3b, "04_output/m3b_lme_varident_collapsed.rds")
saveRDS(m4, "04_output/m4_gamma_glmm.rds")
saveRDS(m4b, "04_output/m4b_gamma_glmm_collapsed.rds")
saveRDS(m5, "04_output/m5_lme_varident_dropped_adhoc.rds")
saveRDS(m5_homog, "04_output/m5_homog_lme.rds")

# Diagnostics
saveRDS(diag_df3, "04_output/m3_diag_df.rds")
saveRDS(re_df3, "04_output/m3_re_df.rds")
saveRDS(diag_df4b, "04_output/m4b_diag_df.rds")
saveRDS(diag_df4b_dev, "04_output/m4b_diag_df_dev.rds")
saveRDS(re_df4b, "04_output/m4b_re_df.rds")

# Comparison table
saveRDS(comparison, "04_output/comparison_table.rds")

# Imputation inputs
saveRDS(top6_beta, "04_output/top6_beta.rds")
saveRDS(pair_map, "04_output/pair_map.rds")

# Working dataset (with collapsed locations, imputed columns, dropped ad hoc)
saveRDS(d, "04_output/d_model_ready.rds")

# Diagnostic plots
ggsave("04_output/m3_diagnostic_plots.png",
       (p1 + p2) / (p3 + p4) + plot_annotation(title = "Model m3: Residual Diagnostics"),
       width = 12, height = 8, dpi = 300)

# LRT result
lrt_result <- anova(m5_homog, m5)
saveRDS(lrt_result, "04_output/lrt_homog_vs_hetero.rds")