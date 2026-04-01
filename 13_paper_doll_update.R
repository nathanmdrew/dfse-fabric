######################################################
### Title: 13_paper_doll_update.R
### Description: Use existing plot logic to recreate before and after imputation
###              paper dolls for the top 4 PAHs, by condition and timing
### Author: Nathan M. Drew (vom8)
### Date: 2026-03-31
#######################################################

# --- Step 0: Setup & Data Prep ------------------------------------------------

library(tidyverse)
library(patchwork)
library(cowplot)

# --- Load saved artifacts from Script 12 --------------------------------------
step1 <- readRDS("12_output/step1_artifacts.rds")
step2 <- readRDS("12_output/step2_artifacts.rds")

top4_cols <- step1$top4_pah_cols

# --- Geometric mean helper (exclude zeros; all-zero → 0) ---------------------
geom_mean <- function(x, na.rm = TRUE) {
  x_pos <- x[x > 0 & !is.na(x)]
  if (length(x_pos) == 0) return(0)
  exp(mean(log(x_pos), na.rm = na.rm))
}

# --- Create as-measured top-4 sum from raw PAH columns ------------------------
d_paper <- step2$d_imp |>
  mutate(
    totalPAH_top4_measured = rowSums(pick(all_of(top4_cols)), na.rm = TRUE)
  ) |>
  select(ParticipantID, Condition, Timing_new, SampleLocation,
         totalPAH_top4_measured, totalPAH_imputed)

# --- QC checks ----------------------------------------------------------------
cat("Rows:", nrow(d_paper), "\n")
cat("Participants:", n_distinct(d_paper$ParticipantID), "\n")
cat("Conditions:", paste(levels(factor(d_paper$Condition)), collapse = ", "), "\n")
cat("Timings:", paste(levels(d_paper$Timing_new), collapse = ", "), "\n")
cat("Locations:", paste(levels(d_paper$SampleLocation), collapse = ", "), "\n")
cat("\nMeasured top-4 summary:\n")
summary(d_paper$totalPAH_top4_measured)
cat("\nImputed top-4 summary:\n")
summary(d_paper$totalPAH_imputed)
cat("\nMeasured zeros:", sum(d_paper$totalPAH_top4_measured == 0), "of", nrow(d_paper),
    "(", round(100 * mean(d_paper$totalPAH_top4_measured == 0), 1), "% )\n")
cat("Imputed zeros:", sum(d_paper$totalPAH_imputed == 0), "of", nrow(d_paper), "\n")

# =============================================================================
# --- Step 1: Body Region Geometry & Core Plotting Function --------------------
# =============================================================================

# --- Body region tile definitions (from 00g_paper_doll_function.R) ------------

body_regions_front <- tribble(
  ~region_name, ~x, ~y, ~width, ~height,
  "Neck",        0,    9,    0.8,  0.6,
  "Chest",       0,    6.75, 2.5,  3,
  "Lower Chest", 0,    4.5,  2.5,  1.5,
  "Sleeve",     -2.5,  5,    1.2,  3,
  "Palm",       -2.5,  2.5,  0.8,  1,
  "Finger",     -2.5,  1.5,  0.8,  0.8,
  "Thumb",      -3.2,  2.5,  0.5,  0.6,
  "Pant",       -0.8,  2,    1.2,  2,
  "Lower Pant", -0.8,  0,    1.2,  1.5,
  "Sock",       -0.8, -1,    1.2,  0.8,
  "Sleeve",      2.5,  5,    1.2,  3,
  "Palm",        2.5,  2.5,  0.8,  1,
  "Finger",      2.5,  1.5,  0.8,  0.8,
  "Thumb",       3.2,  2.5,  0.5,  0.6,
  "Pant",        0.8,  2,    1.2,  2,
  "Lower Pant",  0.8,  0,    1.2,  1.5,
  "Sock",        0.8, -1,    1.2,  0.8,
  "Fly",         0,    3.25, 1.5,  0.5
)

body_regions_back <- tribble(
  ~region_name, ~x, ~y, ~width, ~height,
  "Neck",        0,    9,    0.8,  0.6,
  "Back",        0,    6,    2.5,  4.5,
  "Sleeve",     -2.5,  5,    1.2,  3,
  "Palm",       -2.5,  2.5,  0.8,  1,
  "Finger",     -2.5,  1.5,  0.8,  0.8,
  "Thumb",      -3.2,  2.5,  0.5,  0.6,
  "Pant",       -0.8,  2,    1.2,  2,
  "Lower Pant", -0.8,  0,    1.2,  1.5,
  "Sock",       -0.8, -1,    1.2,  0.8,
  "Sleeve",      2.5,  5,    1.2,  3,
  "Palm",        2.5,  2.5,  0.8,  1,
  "Finger",      2.5,  1.5,  0.8,  0.8,
  "Thumb",       3.2,  2.5,  0.5,  0.6,
  "Pant",        0.8,  2,    1.2,  2,
  "Lower Pant",  0.8,  0,    1.2,  1.5,
  "Sock",        0.8, -1,    1.2,  0.8,
  "Fly",         0,    3.25, 1.5,  0.5
)

# --- Core plotting function (v3) ----------------------------------------------
# Returns a front + back patchwork WITHOUT title or legend.
# Title and legend are handled at the page level.

build_paper_doll <- function(pah_by_location, value_col, max_val,
                             scale_label = "Total PAH\n(µg/g)") {

  shared_scale <- scale_fill_viridis_c(
    option = "plasma",
    limits = c(0, max_val),
    na.value = "grey90",
    name = scale_label,
    begin = 0.05,
    end = 0.95
  )

  doll_theme <- theme_void() +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )

  # --- Front ---
  front_data <- body_regions_front |>
    left_join(pah_by_location, by = c("region_name" = "SampleLocation")) |>
    mutate(fill_val = replace_na(.data[[value_col]], 0),
           fill_val = if_else(fill_val == 0, NA_real_, fill_val))

  p_front <- ggplot(front_data, aes(x = x, y = y, fill = fill_val)) +
    geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
    geom_text(aes(label = region_name), size = 2.5, lineheight = 0.8) +
    shared_scale + coord_fixed() + doll_theme +
    labs(subtitle = "Front")

  # --- Back ---
  back_data <- body_regions_back |>
    left_join(pah_by_location, by = c("region_name" = "SampleLocation")) |>
    mutate(
      fill_val = if_else(region_name == "Back",
                         replace_na(.data[[value_col]], 0), NA_real_),
      fill_val = if_else(fill_val == 0, NA_real_, fill_val)
    )

  p_back <- ggplot(back_data, aes(x = x, y = y, fill = fill_val)) +
    geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
    geom_text(aes(label = if_else(region_name == "Back", "Back", "")),
              size = 2.5, lineheight = 0.8) +
    shared_scale + coord_fixed() + doll_theme +
    labs(subtitle = "Back")

  p_front + p_back
}

# --- Page-level assembly function (v6) ----------------------------------------
build_doll_page <- function(pah_by_location, max_val, page_title) {

  shared_scale <- scale_fill_viridis_c(
    option = "plasma", limits = c(0, max_val), na.value = "grey90",
    name = "Total PAH (µg/g)", begin = 0.05, end = 0.95
  )

  doll_theme <- theme_void() +
    theme(
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "none"
    )

  make_front <- function(val_col) {
    fd <- body_regions_front |>
      left_join(pah_by_location, by = c("region_name" = "SampleLocation")) |>
      mutate(fill_val = replace_na(.data[[val_col]], 0),
             fill_val = if_else(fill_val == 0, NA_real_, fill_val))
    ggplot(fd, aes(x = x, y = y, fill = fill_val)) +
      geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
      geom_text(aes(label = region_name), size = 2.5) +
      shared_scale + coord_fixed() + doll_theme +
      labs(subtitle = "Front")
  }

  make_back <- function(val_col) {
    bd <- body_regions_back |>
      left_join(pah_by_location, by = c("region_name" = "SampleLocation")) |>
      mutate(
        fill_val = if_else(region_name == "Back",
                           replace_na(.data[[val_col]], 0), NA_real_),
        fill_val = if_else(fill_val == 0, NA_real_, fill_val)
      )
    ggplot(bd, aes(x = x, y = y, fill = fill_val)) +
      geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
      geom_text(aes(label = if_else(region_name == "Back", "Back", "")), size = 2.5) +
      shared_scale + coord_fixed() + doll_theme +
      labs(subtitle = "Back")
  }

  # --- Row labels ---
  label_measured <- wrap_elements(
    grid::textGrob("As-Measured", rot = 90,
                   gp = grid::gpar(fontsize = 12, fontface = "bold"))
  )
  label_imputed <- wrap_elements(
    grid::textGrob("Imputed", rot = 90,
                   gp = grid::gpar(fontsize = 12, fontface = "bold"))
  )

  # --- Horizontal separator ---
  h_sep <- wrap_elements(
    grid::linesGrob(x = unit(c(0, 1), "npc"), y = unit(c(0.5, 0.5), "npc"),
                    gp = grid::gpar(col = "grey50", lwd = 1))
  )

  # Layout:
  # A = measured label, B = measured front, C = measured back
  # D = separator (full width)
  # E = imputed label, F = imputed front, G = imputed back
  # H = guide area (full width)
  layout <- "
    ABC
    DDD
    EFG
    HHH
  "

  label_measured + make_front("totalPAH_top4_measured") + make_back("totalPAH_top4_measured") +
    h_sep +
    label_imputed + make_front("totalPAH_imputed") + make_back("totalPAH_imputed") +
    guide_area() +
    plot_layout(
      design = layout,
      widths = c(0.5, 5, 5),
      heights = c(10, 0.3, 10, 1),
      guides = "collect"
    ) +
    plot_annotation(
      title = page_title,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom", legend.key.width = unit(2.5, "cm"))
}

# --- QC test ------------------------------------------------------------------
build_doll_page(test_data, test_max, "AA05 | SS | Doffing")

# Check max across both columns
tibble(
  measured_max = max(d_paper$totalPAH_top4_measured, na.rm = TRUE),
  imputed_max  = max(d_paper$totalPAH_imputed, na.rm = TRUE),
  test_max_used = test_max
)

# =============================================================================
# --- Step 2: Global Color Scale & Data Collapse --------------------------------
# =============================================================================

# --- Global max for consistent color scale ------------------------------------
global_max <- max(d_paper$totalPAH_top4_measured, d_paper$totalPAH_imputed, na.rm = TRUE)
cat("Global max for color scale:", round(global_max, 2), "µg/g\n")

# --- Individual-level: collapse to geomean per participant/condition/timing/location ---
d_individual <- d_paper |>
  summarise(
    totalPAH_top4_measured = geom_mean(totalPAH_top4_measured),
    totalPAH_imputed = geom_mean(totalPAH_imputed),
    .by = c(ParticipantID, Condition, Timing_new, SampleLocation)
  )

cat("\nIndividual dataset:")
cat("\n  Rows:", nrow(d_individual), "(from", nrow(d_paper), ")")
cat("\n  Unique participant/condition/timing combos:",
    n_distinct(paste(d_individual$ParticipantID, d_individual$Condition, d_individual$Timing_new)))

# --- Aggregate-level: geomean across participants per condition/timing/location ---
d_aggregate <- d_paper |>
  summarise(
    totalPAH_top4_measured = geom_mean(totalPAH_top4_measured),
    totalPAH_imputed = geom_mean(totalPAH_imputed),
    .by = c(Condition, Timing_new, SampleLocation)
  )

cat("\n\nAggregate dataset:")
cat("\n  Rows:", nrow(d_aggregate))
cat("\n  Unique condition/timing combos:",
    n_distinct(paste(d_aggregate$Condition, d_aggregate$Timing_new)))

# --- QC: check aggregate ranges vs global max ---------------------------------
cat("\n\nAggregate max (measured):", round(max(d_aggregate$totalPAH_top4_measured), 3))
cat("\nAggregate max (imputed):", round(max(d_aggregate$totalPAH_imputed), 3))
cat("\nGlobal max:", round(global_max, 2))
cat("\n(Aggregate values will use a much smaller portion of the color scale)\n")

# --- Step 2 QC: Test with a different participant to verify -------------------
# Pick someone with more detects for a visually interesting test
detect_counts <- d_individual |>
  summarise(n_detects = sum(totalPAH_top4_measured > 0), .by = ParticipantID) |>
  arrange(desc(n_detects))

cat("\nTop 5 participants by number of detected locations:\n")
print(detect_counts |> head(5))

# Test with the top participant
test_id <- detect_counts$ParticipantID[1]
test_individual <- d_individual |>
  filter(ParticipantID == test_id, Condition == "SS", Timing_new == "Doffing")

build_doll_page(test_individual, global_max,
                paste0(test_id, " | SS | Doffing (Individual)"))

# --- Step 2 QC: Aggregate test plot -------------------------------------------
test_agg <- d_aggregate |>
  filter(Condition == "SS", Timing_new == "Doffing")

# Preview with global scale
build_doll_page(test_agg, global_max, "Aggregate | SS | Doffing (Global Scale)")

# Preview with aggregate-specific scale
agg_max <- max(d_aggregate$totalPAH_top4_measured, d_aggregate$totalPAH_imputed)
build_doll_page(test_agg, agg_max, "Aggregate | SS | Doffing (Aggregate Scale)")

# Check: do all 6 condition/timing combos have all 12 locations?
d_aggregate |>
  count(Condition, Timing_new) |>
  print()

# Which locations are missing?
all_locations <- levels(d_paper$SampleLocation)
d_aggregate |>
  summarise(locations = list(unique(SampleLocation)), .by = c(Condition, Timing_new)) |>
  mutate(
    n_loc = map_int(locations, length),
    missing = map_chr(locations, \(x) paste(setdiff(all_locations, x), collapse = ", "))
  ) |>
  select(-locations)



# --- Page-level assembly function (v7 — handles not-sampled) ------------------
build_doll_page <- function(pah_by_location, max_val, page_title) {

  shared_scale <- scale_fill_viridis_c(
    option = "plasma", limits = c(0, max_val), na.value = "grey90",
    name = "Total PAH (µg/g)", begin = 0.05, end = 0.95
  )

  doll_theme <- theme_void() +
    theme(
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "none"
    )

  make_front <- function(val_col) {
    fd <- body_regions_front |>
      left_join(pah_by_location, by = c("region_name" = "SampleLocation")) |>
      mutate(
        raw_val = .data[[val_col]],
        # Three states: NA = not sampled, 0 = non-detect, >0 = detected
        fill_val = case_when(
          is.na(raw_val) ~ NA_real_,   # not sampled → grey90 via na.value
          raw_val == 0   ~ -Inf,       # non-detect → will map below scale → grey90
          TRUE           ~ raw_val     # detected → color gradient
        ),
        # Use -Inf trick: set to NA so both not-sampled and non-detect get grey90
        fill_val = if_else(fill_val == 0 | is.na(raw_val), NA_real_,
                           if_else(raw_val == 0, NA_real_, raw_val)),
        # Label: region name for sampled locations, "NS" overlay for not-sampled
        not_sampled = is.na(raw_val)
      )

    ggplot(fd, aes(x = x, y = y)) +
      geom_tile(aes(width = width, height = height, fill = fill_val),
                color = "black", linewidth = 0.5) +
      geom_text(aes(label = region_name), size = 2.5,
                data = fd |> filter(!not_sampled)) +
      geom_text(aes(label = paste0(region_name, "\n(NS)")), size = 2, fontface = "italic",
                color = "grey50",
                data = fd |> filter(not_sampled)) +
      shared_scale + coord_fixed() + doll_theme +
      labs(subtitle = "Front")
  }

  make_back <- function(val_col) {
    bd <- body_regions_back |>
      left_join(pah_by_location, by = c("region_name" = "SampleLocation")) |>
      mutate(
        raw_val = .data[[val_col]],
        fill_val = case_when(
          region_name != "Back" ~ NA_real_,   # non-Back regions: no fill
          is.na(raw_val)        ~ NA_real_,   # not sampled
          raw_val == 0           ~ NA_real_,  # non-detect
          TRUE                   ~ raw_val
        ),
        not_sampled = is.na(raw_val) & region_name == "Back",
        is_back = region_name == "Back"
      )

    ggplot(bd, aes(x = x, y = y)) +
      geom_tile(aes(width = width, height = height, fill = fill_val),
                color = "black", linewidth = 0.5) +
      geom_text(aes(label = "Back"), size = 2.5,
                data = bd |> filter(is_back & !not_sampled)) +
      geom_text(aes(label = "Back\n(NS)"), size = 2, fontface = "italic",
                color = "grey50",
                data = bd |> filter(is_back & not_sampled)) +
      shared_scale + coord_fixed() + doll_theme +
      labs(subtitle = "Back")
  }

  # --- Row labels ---
  label_measured <- wrap_elements(
    grid::textGrob("As-Measured", rot = 90,
                   gp = grid::gpar(fontsize = 12, fontface = "bold"))
  )
  label_imputed <- wrap_elements(
    grid::textGrob("Imputed", rot = 90,
                   gp = grid::gpar(fontsize = 12, fontface = "bold"))
  )

  # --- Horizontal separator ---
  h_sep <- wrap_elements(
    grid::linesGrob(x = unit(c(0, 1), "npc"), y = unit(c(0.5, 0.5), "npc"),
                    gp = grid::gpar(col = "grey50", lwd = 1))
  )

  layout <- "
    ABC
    DDD
    EFG
    HHH
  "

  label_measured + make_front("totalPAH_top4_measured") + make_back("totalPAH_top4_measured") +
    h_sep +
    label_imputed + make_front("totalPAH_imputed") + make_back("totalPAH_imputed") +
    guide_area() +
    plot_layout(
      design = layout,
      widths = c(0.5, 5, 5),
      heights = c(10, 0.3, 10, 1),
      guides = "collect"
    ) +
    plot_annotation(
      title = page_title,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom", legend.key.width = unit(2.5, "cm"))
}

# --- Step 2 QC: Test with completed aggregate data ----------------------------
# Use OL/Doffing since it has 2 not-sampled locations (Fly, Lower Chest)
test_agg_ns <- d_aggregate_complete |>
  filter(Condition == "OL", Timing_new == "Doffing")

build_doll_page(test_agg_ns, agg_max, "Aggregate | OL | Doffing")


# --- Complete aggregate data with all locations, distinguishing not-sampled ---
d_aggregate_complete <- d_aggregate |>
  complete(Condition, Timing_new, SampleLocation = all_locations,
           fill = list(totalPAH_top4_measured = NA_real_,
                       totalPAH_imputed = NA_real_))

# Check: should now be 72 rows (6 combos x 12 locations)
cat("Rows:", nrow(d_aggregate_complete), "\n")

# Verify: NA = not sampled, 0 = all non-detect, >0 = has detects
d_aggregate_complete |>
  mutate(status = case_when(
    is.na(totalPAH_top4_measured) ~ "Not Sampled",
    totalPAH_top4_measured == 0   ~ "Non-Detect",
    TRUE                          ~ "Detected"
  )) |>
  count(Condition, Timing_new, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0)

# --- Step 2 QC: Test with completed aggregate data in other conditions/timings ----------------------------

test_agg_ns <- d_aggregate_complete |>
  filter(Condition == "SS", Timing_new == "Doffing")

build_doll_page(test_agg_ns, agg_max, "Aggregate | SS | Doffing")

test_agg_ns <- d_aggregate_complete |>
  filter(Condition == "SL", Timing_new == "Doffing")

build_doll_page(test_agg_ns, agg_max, "Aggregate | SL | Doffing")

test_agg_ns <- d_aggregate_complete |>
  filter(Condition == "OL", Timing_new == "Donning/Firefighting")

build_doll_page(test_agg_ns, agg_max, "Aggregate | OL | Donning/Firefighting")


# =============================================================================
# --- Step 3: Per-Participant Individual PDFs -----------------------------------
# =============================================================================

# Check: do all 23 participants have data for all 3 conditions x 2 timings?
d_individual |>
  distinct(ParticipantID, Condition, Timing_new) |>
  count(ParticipantID, name = "n_combos") |>
  count(n_combos, name = "n_participants")

# If not all 6, which combos are missing?
d_individual |>
  distinct(ParticipantID, Condition, Timing_new) |>
  complete(ParticipantID, Condition, Timing_new) |>
  anti_join(
    d_individual |> distinct(ParticipantID, Condition, Timing_new),
    by = c("ParticipantID", "Condition", "Timing_new")
  )

# --- Identify valid participant/condition/timing combos -----------------------
valid_combos <- d_individual |>
  distinct(ParticipantID, Condition, Timing_new)

cat("Valid combos:", nrow(valid_combos), "of", 23 * 3 * 2, "possible\n")
valid_combos |> count(Condition, name = "n_pages")

# --- Complete individual data with all locations (for not-sampled logic) ------
d_individual_complete <- d_individual |>
  complete(
    nesting(ParticipantID, Condition, Timing_new),
    SampleLocation = all_locations,
    fill = list(totalPAH_top4_measured = NA_real_,
                totalPAH_imputed = NA_real_)
  )

cat("\nIndividual complete rows:", nrow(d_individual_complete), "\n")
cat("Expected:", nrow(valid_combos), "combos x 12 locations =",
    nrow(valid_combos) * 12, "\n")

# --- Create output directory --------------------------------------------------
dir.create("13_output", showWarnings = FALSE)

# --- Generate one PDF per condition -------------------------------------------
conditions <- sort(unique(valid_combos$Condition))

for (cond in conditions) {

  combos_cond <- valid_combos |>
    filter(Condition == cond) |>
    arrange(ParticipantID, Timing_new)

  filename <- paste0("13_output/paper_doll_individual_", cond, ".pdf")
  cat("Generating:", filename, "(", nrow(combos_cond), "pages ) ...\n")

  pdf(filename, width = 10, height = 12)

  for (i in seq_len(nrow(combos_cond))) {

    pid <- combos_cond$ParticipantID[i]
    tmg <- combos_cond$Timing_new[i]

    plot_data <- d_individual_complete |>
      filter(ParticipantID == pid, Condition == cond, Timing_new == tmg)

    page_title <- paste0(pid, " | ", cond, " | ", tmg)

    p <- build_doll_page(plot_data, global_max, page_title)
    print(p)
  }

  dev.off()
  cat("  Saved:", filename, "\n")
}


# =============================================================================
# --- Step 4: Aggregate PDFs ---------------------------------------------------
# =============================================================================

# --- Valid aggregate combos ---------------------------------------------------
valid_agg_combos <- d_aggregate_complete |>
  distinct(Condition, Timing_new) |>
  arrange(Condition, Timing_new)

cat("Aggregate combos:", nrow(valid_agg_combos), "\n")

# --- Single PDF with all 6 pages ---------------------------------------------
filename_agg <- "13_output/paper_doll_aggregate.pdf"
cat("Generating:", filename_agg, "...\n")

pdf(filename_agg, width = 10, height = 12)

for (i in seq_len(nrow(valid_agg_combos))) {

  cond_i <- valid_agg_combos$Condition[i]
  tmg_i  <- valid_agg_combos$Timing_new[i]

  plot_data <- d_aggregate_complete |>
    filter(Condition == cond_i, Timing_new == tmg_i)

  page_title <- paste0("Aggregate (Geometric Mean) | ", cond_i, " | ", tmg_i)

  p <- build_doll_page(plot_data, agg_max, page_title)
  print(p)
}

dev.off()
cat("Saved:", filename_agg, "-", nrow(valid_agg_combos), "pages\n")


# =============================================================================
# --- Step 5: Save Artifacts ---------------------------------------------------
# =============================================================================

artifacts <- list(
  # Data
  d_paper = d_paper,
  d_individual = d_individual,
  d_individual_complete = d_individual_complete,
  d_aggregate = d_aggregate,
  d_aggregate_complete = d_aggregate_complete,
  valid_combos = valid_combos,

  # Scale parameters
  global_max = global_max,
  agg_max = agg_max,

  # Geometry
  body_regions_front = body_regions_front,
  body_regions_back = body_regions_back,

  # Helper
  top4_cols = top4_cols,
  all_locations = all_locations
)

saveRDS(artifacts, "13_output/step13_artifacts.rds")
cat("Saved: 13_output/step13_artifacts.rds\n")
cat("Contents:", paste(names(artifacts), collapse = ", "), "\n")