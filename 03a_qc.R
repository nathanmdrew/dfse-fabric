### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-03-05
### Purpose: Investigate the heavy-right tail of the non-zero new_totalPAH distribution

# ---------------------------------
# --- Step 0: Setup             ---
# ---------------------------------

library(tidyverse)
library(lme4)
library(patchwork)

setwd("C:/Users/vom8/dfse-fabric/")

d <- readRDS("data/cleaned_data.RDS")

d <- d |>
  rename(SampleLocation = `Sample Location`) |>
  mutate(
    SampleLocation = factor(SampleLocation),
    Condition      = factor(Condition),
    Timing_new     = factor(Timing_new),
    ParticipantID  = factor(ParticipantID)
  )

# Load m2 and diagnostics
m2      <- readRDS("03_output/m2_imputed_lmer.rds")
diag_df2 <- readRDS("03_output/m2_diag_df.rds")

# -----------------------------------------
# --- Step 1: Identify the outliers     ---
# -----------------------------------------

# Attach residuals and key variables to the data
d_diag <- d |>
  mutate(
    fitted_m2 = fitted(m2),
    resid_m2  = residuals(m2),
    sresid_m2 = residuals(m2, scaled = TRUE)
  )

# Flag observations with scaled residual > 2
outliers <- d_diag |>
  filter(sresid_m2 > 2) |>
  select(ParticipantID, SampleLocation, Condition, Timing_new,
         new_totalPAH, fitted_m2, sresid_m2) |>
  arrange(desc(sresid_m2))

outliers

# -----------------------------------------
# --- Step 1b: Outlier patterns         ---
# -----------------------------------------

# Tallies by factor
outliers |> count(ParticipantID, sort = TRUE)
outliers |> count(SampleLocation, sort = TRUE)
outliers |> count(Condition, sort = TRUE)
outliers |> count(Timing_new, sort = TRUE)

# Cross-tabulation: Location x Condition
outliers |> count(SampleLocation, Condition, sort = TRUE)

# How many of the 23 participants contribute outliers?
n_distinct(outliers$ParticipantID)

# -----------------------------------------------
# --- Step 2: Examine raw PAH values for outliers
# -----------------------------------------------

# Load top6 beta info for PAH column names
top6_beta <- readRDS("03_output/top6_beta.rds")

# Extract the 6 raw PAH columns + imputed columns for outlier rows
outlier_idx <- which(d_diag$sresid_m2 > 2)

imp_cols <- paste0(top6_beta$pah_name, "_imp")

outlier_pah <- d_diag |>
  slice(outlier_idx) |>
  select(ParticipantID, SampleLocation, Condition, Timing_new,
         all_of(top6_beta$pah_name)) |>
  pivot_longer(
    cols = all_of(top6_beta$pah_name),
    names_to = "PAH",
    values_to = "raw_value"
  )

# How many of the individual PAH values are detected (> 0) vs. non-detect?
outlier_pah |>
  summarise(
    n_total   = n(),
    n_detect  = sum(raw_value > 0, na.rm = TRUE),
    pct_detect = 100 * n_detect / n_total,
    .by = PAH
  )

# Per outlier observation: how many of the 6 PAHs were detected?
outlier_pah |>
  summarise(
    n_detect     = sum(raw_value > 0, na.rm = TRUE),
    total_raw    = sum(raw_value, na.rm = TRUE),
    .by = c(ParticipantID, SampleLocation, Condition, Timing_new)
  ) |>
  arrange(desc(total_raw))

# Is there a pattern with Burn?
d_diag |>
  filter(sresid_m2 > 2) |>
  count(ParticipantID, Burn, sort = TRUE)

# For context, what's the overall distribution of Burn?
d |> count(Burn, sort = TRUE)

# -------------------------------------------
# --- Step 3: Influence diagnostics       ---
# -------------------------------------------

# Recreate the response variable that m2 expects
top6_imp <- readRDS("03_output/top6_beta.rds")

# Rebuild imputed PAH columns and totalPAH_imputed
for (i in seq_len(nrow(top6_imp))) {
  pc   <- top6_imp$pah_col[i]
  lc   <- top6_imp$lod_col[i]
  beta <- top6_imp$beta_mle[i]
  nm   <- top6_imp$pah_name[i]

  pah_vals <- d[[pc]]
  lod_vals <- d[[lc]]
  is_nd    <- !is.na(pah_vals) & !is.na(lod_vals) & (pah_vals <= lod_vals)

  imp_name <- paste0(nm, "_imp")
  d[[imp_name]] <- pah_vals
  d[[imp_name]][is_nd] <- beta * lod_vals[is_nd]
}

imp_cols <- paste0(top6_imp$pah_name, "_imp")
d <- d |>
  mutate(
    totalPAH_imputed = rowSums(pick(all_of(imp_cols)), na.rm = TRUE),
    log_totalPAH_imp = log(totalPAH_imputed)
  )

# Cook's distance and leverage
infl <- influence(m2, "ParticipantID")

# Cook's distance per participant
cooks_d <- cooks.distance(infl)

cooks_df <- tibble(
  ParticipantID = names(cooks_d),
  cooks_d       = as.numeric(cooks_d)
) |>
  arrange(desc(cooks_d))

cooks_df

# Common threshold: 4 / n_groups
threshold <- 4 / n_distinct(d$ParticipantID)

cooks_df |>
  filter(cooks_d > threshold) |>
  mutate(threshold = threshold)


# Check if names carried through
names(cooks_d)

# If NULL, pull from the influence object directly
cooks_df <- tibble(
  ParticipantID = levels(d$ParticipantID),
  cooks_d       = as.numeric(cooks_d)
) |>
  arrange(desc(cooks_d))

cooks_df

# Check the order used by the influence object
unique.del <- unique(model.frame(m2)$ParticipantID)
unique.del

# Rebuild with correct alignment
cooks_df <- tibble(
  ParticipantID = as.character(unique.del),
  cooks_d       = as.numeric(cooks_d)
) |>
  arrange(desc(cooks_d))

cooks_df

# What does AA01 look like?
d |>
  filter(ParticipantID == "AA01") |>
  select(ParticipantID, SampleLocation, Condition, Timing_new,
         totalPAH_imputed, log_totalPAH_imp) |>
  arrange(desc(totalPAH_imputed))

# Which sample locations does only AA01 occupy?
d |>
  count(SampleLocation) |>
  filter(n <= 2)

# -----------------------------------------
# --- Step 4: Visualize outlier patterns ---
# -----------------------------------------

# Residuals colored by key factors

# Residuals vs. fitted, colored by SampleLocation
p1 <- ggplot(d_diag, aes(x = fitted_m2, y = sresid_m2, color = SampleLocation)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(title = "Residuals by Sample Location",
       x = "Fitted values", y = "Scaled residuals") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# Residuals vs. fitted, colored by Condition
p2 <- ggplot(d_diag, aes(x = fitted_m2, y = sresid_m2, color = Condition)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(title = "Residuals by Condition",
       x = "Fitted values", y = "Scaled residuals") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# Residuals vs. fitted, colored by Timing
p3 <- ggplot(d_diag, aes(x = fitted_m2, y = sresid_m2, color = Timing_new)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(title = "Residuals by Timing",
       x = "Fitted values", y = "Scaled residuals") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# Boxplot of residuals by SampleLocation
p4 <- ggplot(d_diag, aes(x = reorder(SampleLocation, sresid_m2, FUN = median),
                          y = sresid_m2)) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Residuals by Sample Location",
       x = NULL, y = "Scaled residuals") +
  theme_minimal(base_size = 12)

(p1 + p2) / (p3 + p4) +
  plot_annotation(title = "Model m2: Outlier Pattern Investigation")

# -----------------------------------------
# --- Step 5: Save artifacts to 03a_output
# -----------------------------------------

# Model diagnostics
saveRDS(d_diag, "03a_output/d_diag.rds")
saveRDS(outliers, "03a_output/outliers.rds")
saveRDS(outlier_pah, "03a_output/outlier_pah.rds")

# Influence diagnostics
saveRDS(infl, "03a_output/influence_m2.rds")
saveRDS(cooks_df, "03a_output/cooks_df.rds")

# Diagnostic plots
ggsave("03a_output/m2_outlier_patterns.png",
       (p1 + p2) / (p3 + p4) + plot_annotation(title = "Model m2: Outlier Pattern Investigation"),
       width = 14, height = 8, dpi = 300)