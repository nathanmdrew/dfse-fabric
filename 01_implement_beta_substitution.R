### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-03-13
### Purpose: Begin implementing beta-substitution

library(tidyverse)
library(tibble)

setwd("C:/Users/vom8/dfse-fabric/")

d <- readRDS(file="data/cleaned_data.RDS")

### Quick EDA for LODs
# 8:22 are PAH
# 24:38 are LODs

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

# Row-level exact equality check
eq_rows <- map_dfr(seq_len(nrow(pair_map)), function(i) {
  pcol <- pair_map$pah_col[i]
  lcol <- pair_map$lod_col[i]
  
  out <- d %>%
    transmute(
      row_id = row_number(),
      pair_id = pair_map$pair_id[i],
      pah_name = pair_map$pah_name[i],
      lod_name = pair_map$lod_name[i],
      pah_value = .[[pcol]],
      lod_value = .[[lcol]],
      equal_exact = !is.na(pah_value) & !is.na(lod_value) & (pah_value == lod_value)
    ) %>%
    filter(equal_exact)
  
  out
})

# Summary by PAH
eq_summary <- eq_rows %>%
  count(pair_id, pah_name, lod_name, name = "n_equal") %>%
  arrange(desc(n_equal))

eq_summary
eq_rows

# !!! So there are 10 records with measurements = LOD
# !!! Decision rule - if value <= LOD, substitute.
# !!! Could do sensitivity or other comparison to value < LOD instead

# Rank PAHs, descending by non-missingness/not below LoD
pah_detect_rank <- pair_map %>%
  mutate(
    n_total = nrow(d),
    n_detect = map2_int(
      pah_col, lod_col,
      ~sum(!is.na(d[[.x]]) & !is.na(d[[.y]]) & d[[.x]] > d[[.y]])
    ),
    n_nondetect = map2_int(
      pah_col, lod_col,
      ~sum(!is.na(d[[.x]]) & !is.na(d[[.y]]) & d[[.x]] <= d[[.y]])
    ),
    n_zero = map_int(pah_col, ~sum(!is.na(d[[.x]]) & d[[.x]] == 0)),
    pct_detect = 100 * n_detect / n_total
  ) %>%
  select(pair_id, pah_name, lod_name, n_total, n_detect, n_nondetect, n_zero, pct_detect) %>%
  arrange(desc(n_detect), desc(pct_detect), pah_name)

pah_detect_rank

# Top 6 here (NB: not only after doffing like the paper) still match the paper
# 5 PAHs are totally not detected.

names(d)

# 1) Set these to your actual column names in `d`
site_col   <- "Sample Location"      # e.g., "sampling_site", "sample_site", "location"
period_col <- "Timing"  # e.g., "collection_period", "timing", "phase"

# 2) Normalize text and subset to glove locations + doffing
d_glove_doff <- d %>%
  mutate(
    site_norm   = str_to_lower(as.character(.data[[site_col]])),
    period_norm = str_to_lower(as.character(.data[[period_col]]))
  ) %>%
  filter(str_detect(site_norm, "palm|finger|thumb")) %>%
  filter(str_detect(period_norm, "doff")) %>%
  select(-site_norm, -period_norm)

# Optional checks
d %>% count(.data[[site_col]], .data[[period_col]], sort = TRUE)
d_glove_doff %>% count(.data[[site_col]], .data[[period_col]], sort = TRUE)

#------------------- Example beta substitution for one PAH -------------------
# Beta-substitution for a single PAH
# Reference: Ganser & Hewett - beta is estimated from detected observations

beta_substitute_pah <- function(data, pah_col_idx, lod_col_idx, pah_label = NULL) {
  
  pah_vals <- data[[pah_col_idx]]
  lod_vals <- data[[lod_col_idx]]
  
  # Classify rows
  has_both   <- !is.na(pah_vals) & !is.na(lod_vals)
  is_detect  <- has_both & pah_vals > lod_vals
  is_nd      <- has_both & pah_vals <= lod_vals
  
  # Detected values only
  x_det <- pah_vals[is_detect]
  l_det <- lod_vals[is_detect]
  
  # Estimate beta: mean ratio of detected value to its paired LOD
  # (a simple moment-based estimate; Ganser & Hewett use MLE of lognormal
  #  truncated at LOD -- this is the approximation form)
  beta_hat <- mean(x_det / l_det)
  
  # Apply substitution: imputed value = beta_hat * LOD for non-detects
  imputed <- ifelse(is_nd, beta_hat * lod_vals, pah_vals)
  
  tibble(
    pah_label   = pah_label,
    pah_raw     = pah_vals,
    lod_val     = lod_vals,
    is_detect   = is_detect,
    is_nd       = is_nd,
    beta_hat    = beta_hat,
    imputed_val = imputed
  )
}

# Run for Dibenzo(a,h)anthracene (pair_id = 9, pah_col = 17, lod_col = 31)
target <- pair_map |> filter(pair_id == 9)

result_dibenzo <- beta_substitute_pah(
  data        = d,
  pah_col_idx = target$pah_col,
  lod_col_idx = target$lod_col,
  pah_label   = target$pah_name
)

result_dibenzo

# Build three-group comparison: Detected (raw), Non-detect (β-imputed), LOD
lod_dat <- result_dibenzo |>
  filter(!is.na(lod_val)) |>
  transmute(
    pah_label = pah_label,
    value     = lod_val,
    group     = "LOD",
    is_detect = FALSE,
    is_nd     = FALSE,
    beta_hat  = NA_real_
  )

# Prepare plot data: detected raw values vs imputed ND values
plot_dat <- result_dibenzo |>
  filter(!is.na(pah_raw) & !is.na(lod_val)) |>
  mutate(
    value = if_else(is_detect, pah_raw, imputed_val),
    group = case_when(
      is_detect ~ "Detected (raw)",
      is_nd     ~ "Non-detect (β-imputed)",
      TRUE      ~ NA_character_
    )
  ) |>
  filter(!is.na(group))

plot_combined <- plot_dat |>
  select(pah_label, value, group, is_detect, is_nd, beta_hat) |>
  bind_rows(lod_dat)

ggplot(plot_combined, aes(x = value, fill = group, color = group)) +
  geom_density(alpha = 0.30, linewidth = 0.8) +
  geom_rug(
    data = filter(plot_combined, group != "LOD"),
    alpha = 0.4, length = unit(0.03, "npc")
  ) +
  scale_fill_manual(values = c(
    "Detected (raw)"         = "#E69F00",
    "Non-detect (β-imputed)" = "#0072B2",
    "LOD"                    = "#CC79A7"
  )) +
  scale_color_manual(values = c(
    "Detected (raw)"         = "#E69F00",
    "Non-detect (β-imputed)" = "#0072B2",
    "LOD"                    = "#CC79A7"
  )) +
  labs(
    title    = "Dibenzo(a,h)anthracene: Detected, β-Imputed, and LOD Distributions",
    subtitle = paste0(
      "β̂ = ", round(unique(plot_dat$beta_hat), 2),
      " | n detected = ", sum(plot_dat$is_detect),
      " | n imputed = ",  sum(plot_dat$is_nd)
    ),
    x     = "Concentration",
    y     = "Density",
    fill  = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

#----------------------- Look for LOD clusters -----------------------
result_dibenzo |>
  filter(!is.na(lod_val)) |>
  count(lod_val, sort = TRUE) |>
  mutate(pct = 100 * n / sum(n))

result_dibenzo |>
  filter(!is.na(lod_val)) |>
  ggplot(aes(x = lod_val)) +
  geom_histogram(binwidth = 0.05, fill = "#CC79A7", color = "white", alpha = 0.8) +
  geom_rug(alpha = 0.4, length = unit(0.03, "npc")) +
  labs(
    title = "LOD Distribution: Dibenzo(a,h)anthracene",
    x     = "LOD Value",
    y     = "Count"
  ) +
  theme_minimal(base_size = 13)

#----------------------- Implement the MLE approach for beta estimation -----------------------

library(EnvStats)

# --- 1. Extract observed values and censoring indicator from result_dibenzo ---
# Rows with both PAH and LOD values present
mle_dat <- result_dibenzo |>
  filter(!is.na(pah_raw) & !is.na(lod_val))

obs_vals_real <- mle_dat$pah_raw   # detected = true value; ND = LOD value (as recorded)
is_nd_real    <- mle_dat$is_nd     # TRUE = left-censored at LOD

cat("n total (non-missing):", nrow(mle_dat), "\n")
cat("n detected:           ", sum(!is_nd_real), "\n")
cat("n non-detect:         ", sum(is_nd_real), "\n")
cat("Empirical censoring:  ", round(100 * mean(is_nd_real), 1), "%\n")

# Check data for issues
mle_dat |>
  filter(is_nd) |>
  count(pah_raw) |>
  arrange(pah_raw)

obs_vals_real <- if_else(mle_dat$is_nd, mle_dat$lod_val, mle_dat$pah_raw)
is_nd_real    <- mle_dat$is_nd

# --- 2. Fit censored lognormal via MLE ---
fit_real <- elnormCensored(
  x        = obs_vals_real,
  censored = is_nd_real,
  method   = "mle"
)

# First diagnostic check
mu_hat_real    <- fit_real$parameters[["meanlog"]]
sigma_hat_real <- fit_real$parameters[["sdlog"]]

cat("MLE mu_hat:", round(mu_hat_real, 3), "| sigma_hat:", round(sigma_hat_real, 3), "\n")

# --- 3. Derive MLE beta using median LOD as representative censoring point ---
# Beta = E[X | X <= lod_ref] / lod_ref  (truncated lognormal expectation)
lod_ref_real <- median(mle_dat$lod_val)

z_real        <- (log(lod_ref_real) - mu_hat_real) / sigma_hat_real
beta_hat_mle  <- exp(mu_hat_real + sigma_hat_real^2 / 2) *
                   pnorm(z_real - sigma_hat_real) / pnorm(z_real) / lod_ref_real

cat("MLE beta_hat:", round(beta_hat_mle, 3), "\n")
cat("Naive beta_hat (from result_dibenzo):", round(unique(result_dibenzo$beta_hat[!is.na(result_dibenzo$beta_hat)]), 3), "\n")

# --- 4. Apply MLE beta substitution ---
mle_dat <- mle_dat |>
  mutate(
    imputed_mle = if_else(is_nd, beta_hat_mle * lod_val, pah_raw)
  )

# --- 5. Compare naive vs. MLE imputation side-by-side ---
# Summary statistics: GM, GSD, P95 for detected-only, naive-imputed, MLE-imputed
gm_fn  <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))
gsd_fn <- function(x) exp(sd(log(x[x > 0]),   na.rm = TRUE))

compare_metrics <- tribble(
  ~dataset,          ~values,
  "Detected only",   mle_dat$pah_raw[!is_nd_real],
  "Naive β-imputed", if_else(mle_dat$is_nd, mle_dat$imputed_val, mle_dat$pah_raw),
  "MLE β-imputed",   mle_dat$imputed_mle
) |>
  mutate(
    GM  = map_dbl(values, gm_fn),
    GSD = map_dbl(values, gsd_fn),
    P95 = map_dbl(values, ~quantile(.x, 0.95, na.rm = TRUE)),
    n   = map_int(values, ~sum(!is.na(.x)))
  ) |>
  select(-values)

compare_metrics

# --- 6. KDE overlay: Detected, Naive β, MLE β ---
kde_plot_dat <- bind_rows(
  mle_dat |> filter(!is_nd) |>
    transmute(value = pah_raw,     distribution = "Detected (raw)"),
  mle_dat |>
    transmute(value = if_else(is_nd, imputed_val, pah_raw), distribution = "Naive β-imputed"),
  mle_dat |>
    transmute(value = imputed_mle,  distribution = "MLE β-imputed")
)

ggplot(kde_plot_dat, aes(x = value, fill = distribution, color = distribution)) +
  geom_density(alpha = 0.30, linewidth = 0.8) +
  scale_fill_manual(values = c(
    "Detected (raw)"   = "#E69F00",
    "Naive β-imputed"  = "#0072B2",
    "MLE β-imputed"    = "#009E73"
  )) +
  scale_color_manual(values = c(
    "Detected (raw)"   = "#E69F00",
    "Naive β-imputed"  = "#0072B2",
    "MLE β-imputed"    = "#009E73"
  )) +
  labs(
    title    = "Dibenzo(a,h)anthracene: Naive vs. MLE β-Substitution",
    subtitle = paste0(
      "Naive β̂ = ", round(unique(result_dibenzo$beta_hat[!is.na(result_dibenzo$beta_hat)]), 3),
      " | MLE β̂ = ", round(beta_hat_mle, 3),
      " | Censoring = ", round(100 * mean(is_nd_real), 1), "%",
      " | n = ", nrow(mle_dat)
    ),
    x = "Concentration", y = "Density",
    fill = NULL, color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

#---------------- ECDF plots----------------
ecdf_dat <- bind_rows(
  mle_dat |> transmute(value = if_else(is_nd, imputed_val, pah_raw), distribution = "Naive β-imputed"),
  mle_dat |> transmute(value = imputed_mle,                           distribution = "MLE β-imputed"),
  mle_dat |> filter(!is_nd) |> transmute(value = pah_raw,            distribution = "Detected (raw)")
)

ggplot(ecdf_dat, aes(x = value, color = distribution)) +
  stat_ecdf(linewidth = 0.9) +
  scale_color_manual(values = c(
    "Detected (raw)"   = "#E69F00",
    "Naive β-imputed"  = "#0072B2",
    "MLE β-imputed"    = "#009E73"
  )) +
  labs(
    title = "Dibenzo(a,h)anthracene: ECDF Comparison",
    x = "Concentration", y = "Cumulative Probability",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

#----------------------- ECDF plots with the as-is line
ecdf_dat2 <- bind_rows(
  mle_dat |> transmute(value = pah_raw,                                          distribution = "Raw data (incl. zeros)"),
  mle_dat |> transmute(value = if_else(is_nd, imputed_val, pah_raw),             distribution = "Naive β-imputed"),
  mle_dat |> transmute(value = imputed_mle,                                      distribution = "MLE β-imputed"),
  mle_dat |> filter(!is_nd) |> transmute(value = pah_raw,                        distribution = "Detected (raw)")
)

ggplot(ecdf_dat2, aes(x = value, color = distribution)) +
  stat_ecdf(linewidth = 0.9) +
  scale_color_manual(values = c(
    "Raw data (incl. zeros)" = "#CC79A7",
    "Detected (raw)"         = "#E69F00",
    "Naive β-imputed"        = "#0072B2",
    "MLE β-imputed"          = "#009E73"
  )) +
  labs(
    title    = "Dibenzo(a,h)anthracene: ECDF Comparison",
    subtitle = paste0("n = ", nrow(mle_dat), " | ", sum(mle_dat$is_nd), " non-detects (", round(100 * mean(mle_dat$is_nd), 1), "%) recorded as zero"),
    x = "Concentration", y = "Cumulative Probability",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
