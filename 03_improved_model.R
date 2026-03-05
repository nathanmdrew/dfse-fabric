### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-03-05
### Purpose: Investigate imputing values and refitting the model
### Steps


# Required packages (install once):
#install.packages("tidyverse")
#install.packages("EnvStats")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("emmeans")
#install.packages("patchwork")

library(tidyverse)
library(EnvStats)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)

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

# Fit LMM with same structure as baseline m1
m2 <- lmer(
  log_totalPAH_imp ~ Condition + Timing_new + SampleLocation + (1 | ParticipantID),
  data = d
)

summary(m2)

# ---------------------------------
# --- Step 5: Residual diagnostics
# ---------------------------------

diag_df2 <- tibble(
  fitted = fitted(m2),
  resid  = residuals(m2),
  sresid = residuals(m2, scaled = TRUE)
)

# 1. Residuals vs. fitted
p1 <- ggplot(diag_df2, aes(x = fitted, y = sresid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, color = "steelblue", linewidth = 0.8) +
  labs(title = "Residuals vs. Fitted", x = "Fitted values", y = "Scaled residuals") +
  theme_minimal(base_size = 12)

# 2. QQ plot of residuals
p2 <- ggplot(diag_df2, aes(sample = sresid)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(color = "red") +
  labs(title = "Normal Q-Q", x = "Theoretical quantiles", y = "Scaled residuals") +
  theme_minimal(base_size = 12)

# 3. Histogram of residuals
p3 <- ggplot(diag_df2, aes(x = sresid)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "black") +
  labs(title = "Residual Distribution", x = "Scaled residuals", y = "Count") +
  theme_minimal(base_size = 12)

# 4. Random effects QQ
re_df2 <- tibble(
  re = ranef(m2)$ParticipantID[["(Intercept)"]]
)

p4 <- ggplot(re_df2, aes(sample = re)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(color = "red") +
  labs(title = "Random Effects Q-Q", x = "Theoretical quantiles",
       y = "Participant intercepts") +
  theme_minimal(base_size = 12)

# Combine
(p1 + p2) / (p3 + p4) +
  plot_annotation(title = "Model m2: Residual Diagnostics")

# -------------------------------------------------------
# --- Step 6: Variance components, ANOVA, and contrasts
# -------------------------------------------------------

# Variance components / ICC
vc2 <- as.data.frame(VarCorr(m2))
icc2 <- vc2$vcov[1] / sum(vc2$vcov)

tibble(
  component = c("ParticipantID", "Residual"),
  variance  = vc2$vcov,
  std_dev   = vc2$sdcor,
  ICC       = c(icc2, NA_real_)
)

# Type III F-tests
anova(m2, type = 3)

# Condition contrasts
emm_cond2 <- emmeans(m2, "Condition")
emm_cond2

contrast(emm_cond2, method = "pairwise", adjust = "tukey")

# -----------------------------------------------
# --- Step 7: Save artifacts to 03_output     ---
# -----------------------------------------------

# Model objects
saveRDS(m2, "03_output/m2_imputed_lmer.rds")

# Imputation inputs
saveRDS(pair_map, "03_output/pair_map.rds")
saveRDS(top6_beta, "03_output/top6_beta.rds")

# Comparison table
saveRDS(comparison, "03_output/comparison_table.rds")

# Diagnostics
saveRDS(diag_df2, "03_output/m2_diag_df.rds")
saveRDS(re_df2, "03_output/m2_re_df.rds")

# Variance components and ICC
saveRDS(vc2, "03_output/m2_variance_components.rds")
saveRDS(icc2, "03_output/m2_icc.rds")

# Contrasts
saveRDS(emm_cond2, "03_output/m2_emmeans_condition.rds")

# Diagnostic plots
ggsave("03_output/m2_diagnostic_plots.png",
       (p1 + p2) / (p3 + p4) + plot_annotation(title = "Model m2: Residual Diagnostics"),
       width = 12, height = 8, dpi = 300)

