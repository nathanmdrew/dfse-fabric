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
