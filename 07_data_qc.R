# ============================================================
# 07_data_qc.R
# Investigate multiple measures per Timing_new
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)

setwd("C:/Users/vom8/dfse-fabric/")

d <- readRDS(file="data/cleaned_data.RDS")

# -----------------------------------
# --- Step 0: Setup and data prep ---
# -----------------------------------

# 0a) Confirm required columns exist
required_cols <- c("ParticipantID", "Condition", "Sample Location", "Timing_new", "rownum")
missing_cols <- setdiff(required_cols, names(d))
stopifnot(length(missing_cols) == 0)

# 0b) Build a clean analysis frame for timing-multiplicity checks
d_timing <- d |>
  transmute(
    rownum,
    ParticipantID,
    Condition,
    SampleLocation = `Sample Location`,
    Timing_new,
    SampleID,
    SampleID2,
    new_totalPAH
  ) |>
  mutate(
    across(c(ParticipantID, Condition, SampleLocation, Timing_new), as.character)
  )

# 0c) Basic QC summary (used in Step 1)
step0_qc <- d_timing |>
  summarize(
    n_rows = n(),
    n_participants = n_distinct(ParticipantID),
    n_condition = n_distinct(Condition),
    n_locations = n_distinct(SampleLocation),
    n_missing_timing_new = sum(is.na(Timing_new)),
    n_missing_any_key = sum(
      is.na(ParticipantID) | is.na(Condition) | is.na(SampleLocation) | is.na(Timing_new)
    )
  )

list(
  d_timing = d_timing,
  step0_qc = step0_qc
)


# ---------------------------------------------
# --- Step 1: Quantify multiplicity patterns ---
# ---------------------------------------------

# 1a) Multiplicity in raw row count within ParticipantID/Condition/SampleLocation
step1_rows_by_group <- d_timing |>
  count(ParticipantID, Condition, SampleLocation, name = "n_rows_group") |>
  mutate(has_multiple_rows = n_rows_group > 1L)

step1_rows_summary <- step1_rows_by_group |>
  summarize(
    n_groups = n(),
    n_groups_with_multiple_rows = sum(has_multiple_rows),
    pct_groups_with_multiple_rows = mean(has_multiple_rows),
    total_extra_rows = sum(pmax(n_rows_group - 1L, 0L))
  )

# 1b) Multiplicity in distinct Timing_new values within ParticipantID/Condition/SampleLocation
step1_timing_by_group <- d_timing |>
  distinct(ParticipantID, Condition, SampleLocation, Timing_new) |>
  count(ParticipantID, Condition, SampleLocation, name = "n_distinct_timing") |>
  mutate(has_multiple_timing = n_distinct_timing > 1L)

step1_timing_summary <- step1_timing_by_group |>
  summarize(
    n_groups = n(),
    n_groups_with_multiple_timing = sum(has_multiple_timing),
    pct_groups_with_multiple_timing = mean(has_multiple_timing)
  )

# 1c) Where this happens most often (Condition x SampleLocation)
step1_rows_by_condition_location <- step1_rows_by_group |>
  summarize(
    n_groups = n(),
    n_groups_with_multiple_rows = sum(has_multiple_rows),
    pct_groups_with_multiple_rows = mean(has_multiple_rows),
    total_extra_rows = sum(pmax(n_rows_group - 1L, 0L)),
    .by = c(Condition, SampleLocation)
  ) |>
  arrange(desc(pct_groups_with_multiple_rows), desc(total_extra_rows))

step1_timing_by_condition_location <- step1_timing_by_group |>
  summarize(
    n_groups = n(),
    n_groups_with_multiple_timing = sum(has_multiple_timing),
    pct_groups_with_multiple_timing = mean(has_multiple_timing),
    .by = c(Condition, SampleLocation)
  ) |>
  arrange(desc(pct_groups_with_multiple_timing))

# 1d) Detailed groups to inspect
step1_multiple_rows_detail <- step1_rows_by_group |>
  filter(has_multiple_rows) |>
  arrange(desc(n_rows_group), ParticipantID, Condition, SampleLocation)

step1_multiple_timing_detail <- step1_timing_by_group |>
  filter(has_multiple_timing) |>
  arrange(desc(n_distinct_timing), ParticipantID, Condition, SampleLocation)

list(
  rows_summary = step1_rows_summary,
  timing_summary = step1_timing_summary,
  rows_by_condition_location = step1_rows_by_condition_location,
  timing_by_condition_location = step1_timing_by_condition_location,
  multiple_rows_detail = step1_multiple_rows_detail,
  multiple_timing_detail = step1_multiple_timing_detail
)


# ---------------------------------------------------------------
# --- Step 1 (revised): Duplicates within ID/Cond/Loc/Timing  ---
# ---------------------------------------------------------------

# 1a) Count rows per full key
step1_fullkey_counts <- d_timing |>
  filter(
    !is.na(ParticipantID),
    !is.na(Condition),
    !is.na(SampleLocation),
    !is.na(Timing_new)
  ) |>
  count(
    ParticipantID, Condition, SampleLocation, Timing_new,
    name = "n_rows_fullkey"
  ) |>
  mutate(has_duplicate_fullkey = n_rows_fullkey > 1L)

# 1b) Overall duplicate burden
step1_fullkey_summary <- step1_fullkey_counts |>
  summarize(
    n_fullkey_groups = n(),
    n_groups_with_duplicates = sum(has_duplicate_fullkey),
    pct_groups_with_duplicates = mean(has_duplicate_fullkey),
    total_extra_rows = sum(pmax(n_rows_fullkey - 1L, 0L))
  )

# 1c) Where duplicates occur most (Condition x SampleLocation x Timing_new)
step1_fullkey_by_cls <- step1_fullkey_counts |>
  summarize(
    n_groups = n(),
    n_groups_with_duplicates = sum(has_duplicate_fullkey),
    pct_groups_with_duplicates = mean(has_duplicate_fullkey),
    total_extra_rows = sum(pmax(n_rows_fullkey - 1L, 0L)),
    .by = c(Condition, SampleLocation, Timing_new)
  ) |>
  arrange(desc(pct_groups_with_duplicates), desc(total_extra_rows))

# 1d) Exact duplicated groups for inspection
step1_fullkey_duplicate_detail <- step1_fullkey_counts |>
  filter(has_duplicate_fullkey) |>
  arrange(desc(n_rows_fullkey), ParticipantID, Condition, SampleLocation, Timing_new)

list(
  fullkey_summary = step1_fullkey_summary,
  fullkey_by_condition_location_timing = step1_fullkey_by_cls,
  fullkey_duplicate_detail = step1_fullkey_duplicate_detail
)

qc <- d |> filter(ParticipantID=="BA23" & Condition=="OL" 
                 & `Sample Location`=="Palm" & Timing_new=="Donning/Firefighting")
qc #yeah, 3 measures on Palm.

unique(step1_fullkey_duplicate_detail$ParticipantID) # 12


# ---------------------------------------------
# --- Step 2: Save Step 0/1 artifacts as RDS ---
# ---------------------------------------------


artifact_dir <- "07_output"
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)

# Add/remove names as needed
artifact_names <- c(
  "d_timing",
  "step0_qc",
  "step1_rows_by_group",
  "step1_rows_summary",
  "step1_rows_by_condition_location",
  "step1_multiple_rows_detail",
  "step1_timing_by_group",
  "step1_timing_summary",
  "step1_timing_by_condition_location",
  "step1_multiple_timing_detail",
  "step1_fullkey_counts",
  "step1_fullkey_summary",
  "step1_fullkey_by_cls",
  "step1_fullkey_duplicate_detail"
)

artifact_names_present <- artifact_names[artifact_names %in% ls()]

walk(
  artifact_names_present,
  \(obj_name) {
    saveRDS(
      get(obj_name, inherits = TRUE),
      file = file.path(artifact_dir, paste0(obj_name, ".rds"))
    )
  }
)

saved_artifacts <- tibble::tibble(
  object = artifact_names_present,
  path = file.path(artifact_dir, paste0(artifact_names_present, ".rds"))
)

saved_artifacts