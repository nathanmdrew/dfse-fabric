# ============================================================
# 06_paper_doll_imputed.R
# Paper doll visualizations of MLE-imputed top 6 PAH data
# ============================================================

# ---------------------------------
# --- Step 0: Setup             ---
# ---------------------------------

library(tidyverse)
library(EnvStats)
library(patchwork)
library(cowplot)

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

# Collapse Shirt → Sleeve only (keep Fly and Lower Chest)
d <- d |>
  mutate(
    SampleLocation = fct_recode(SampleLocation, "Sleeve" = "Shirt") |> fct_drop()
  )

# Load top 6 betas and apply MLE imputation
top6_beta <- readRDS("04_output/top6_beta.rds")

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
  mutate(totalPAH_imputed = rowSums(pick(all_of(imp_cols)), na.rm = TRUE))

# Compute global min/max for consistent color scale across all documents
global_min <- min(d$totalPAH_imputed, na.rm = TRUE)
global_max <- max(d$totalPAH_imputed, na.rm = TRUE)

# Verify
d |> count(SampleLocation, sort = TRUE)
d |> count(Condition)
d |> count(Timing_new)
cat("Global range:", round(global_min, 4), "to", round(global_max, 4), "\n")


# -----------------------------------------
# --- Step 1: Define body region geometry --
# -----------------------------------------

# Front view body regions
front_regions <- tribble(
  ~region,        ~x,    ~y,
  # Head
  "Head",         4.5,   16.5,
  "Head",         5.5,   16.5,
  "Head",         5.5,   17.5,
  "Head",         4.5,   17.5,
  # Neck
  "Neck",         4.5,   15.5,
  "Neck",         5.5,   15.5,
  "Neck",         5.5,   16.5,
  "Neck",         4.5,   16.5,
  # Chest
  "Chest",        3.5,   13.0,
  "Chest",        6.5,   13.0,
  "Chest",        6.5,   15.5,
  "Chest",        3.5,   15.5,
  # Lower Chest
  "Lower Chest",  3.5,   11.0,
  "Lower Chest",  6.5,   11.0,
  "Lower Chest",  6.5,   13.0,
  "Lower Chest",  3.5,   13.0,
  # Fly
  "Fly",          4.5,   9.0,
  "Fly",          5.5,   9.0,
  "Fly",          5.5,   11.0,
  "Fly",          4.5,   11.0,
  # Left Sleeve
  "Sleeve",       1.5,   13.0,
  "Sleeve",       3.5,   13.0,
  "Sleeve",       3.5,   15.5,
  "Sleeve",       1.5,   15.5,
  # Right Sleeve
  "Sleeve_R",     6.5,   13.0,
  "Sleeve_R",     8.5,   13.0,
  "Sleeve_R",     8.5,   15.5,
  "Sleeve_R",     6.5,   15.5,
  # Left Pant
  "Pant",         3.5,   5.0,
  "Pant",         4.5,   5.0,
  "Pant",         4.5,   9.0,
  "Pant",         3.5,   9.0,
  # Right Pant
  "Pant_R",       5.5,   5.0,
  "Pant_R",       6.5,   5.0,
  "Pant_R",       6.5,   9.0,
  "Pant_R",       5.5,   9.0,
  # Left Lower Pant
  "Lower Pant",   3.5,   1.0,
  "Lower Pant",   4.5,   1.0,
  "Lower Pant",   4.5,   5.0,
  "Lower Pant",   3.5,   5.0,
  # Right Lower Pant
  "Lower Pant_R", 5.5,   1.0,
  "Lower Pant_R", 6.5,   1.0,
  "Lower Pant_R", 6.5,   5.0,
  "Lower Pant_R", 5.5,   5.0,
  # Left Sock
  "Sock",         3.5,   0.0,
  "Sock",         4.5,   0.0,
  "Sock",         4.5,   1.0,
  "Sock",         3.5,   1.0,
  # Right Sock
  "Sock_R",       5.5,   0.0,
  "Sock_R",       6.5,   0.0,
  "Sock_R",       6.5,   1.0,
  "Sock_R",       5.5,   1.0,
  # Left Palm
  "Palm",         0.5,   11.0,
  "Palm",         1.5,   11.0,
  "Palm",         1.5,   12.5,
  "Palm",         0.5,   12.5,
  # Right Palm
  "Palm_R",       8.5,   11.0,
  "Palm_R",       9.5,   11.0,
  "Palm_R",       9.5,   12.5,
  "Palm_R",       8.5,   12.5,
  # Left Finger
  "Finger",       0.5,   12.5,
  "Finger",       1.5,   12.5,
  "Finger",       1.5,   14.0,
  "Finger",       0.5,   14.0,
  # Right Finger
  "Finger_R",     8.5,   12.5,
  "Finger_R",     9.5,   12.5,
  "Finger_R",     9.5,   14.0,
  "Finger_R",     8.5,   14.0,
  # Left Thumb
  "Thumb",        0.0,   12.0,
  "Thumb",        0.5,   12.0,
  "Thumb",        0.5,   13.0,
  "Thumb",        0.0,   13.0,
  # Right Thumb
  "Thumb_R",      9.5,   12.0,
  "Thumb_R",      10.0,  12.0,
  "Thumb_R",      10.0,  13.0,
  "Thumb_R",      9.5,   13.0
)

# Back view body regions
back_regions <- tribble(
  ~region,        ~x,    ~y,
  # Head
  "Head_B",       4.5,   16.5,
  "Head_B",       5.5,   16.5,
  "Head_B",       5.5,   17.5,
  "Head_B",       4.5,   17.5,
  # Neck
  "Neck_B",       4.5,   15.5,
  "Neck_B",       5.5,   15.5,
  "Neck_B",       5.5,   16.5,
  "Neck_B",       4.5,   16.5,
  # Back
  "Back",         3.5,   11.0,
  "Back",         6.5,   11.0,
  "Back",         6.5,   15.5,
  "Back",         3.5,   15.5,
  # Left Sleeve
  "Sleeve_B",     1.5,   13.0,
  "Sleeve_B",     3.5,   13.0,
  "Sleeve_B",     3.5,   15.5,
  "Sleeve_B",     1.5,   15.5,
  # Right Sleeve
  "Sleeve_BR",    6.5,   13.0,
  "Sleeve_BR",    8.5,   13.0,
  "Sleeve_BR",    8.5,   15.5,
  "Sleeve_BR",    6.5,   15.5,
  # Left Pant
  "Pant_B",       3.5,   5.0,
  "Pant_B",       4.5,   5.0,
  "Pant_B",       4.5,   11.0,
  "Pant_B",       3.5,   11.0,
  # Right Pant
  "Pant_BR",      5.5,   5.0,
  "Pant_BR",      6.5,   5.0,
  "Pant_BR",      6.5,   11.0,
  "Pant_BR",      5.5,   11.0,
  # Left Lower Pant
  "Lower Pant_B", 3.5,   1.0,
  "Lower Pant_B", 4.5,   1.0,
  "Lower Pant_B", 4.5,   5.0,
  "Lower Pant_B", 3.5,   5.0,
  # Right Lower Pant
  "Lower Pant_BR",5.5,   1.0,
  "Lower Pant_BR",6.5,   1.0,
  "Lower Pant_BR",6.5,   5.0,
  "Lower Pant_BR",5.5,   5.0,
  # Left Sock
  "Sock_B",       3.5,   0.0,
  "Sock_B",       4.5,   0.0,
  "Sock_B",       4.5,   1.0,
  "Sock_B",       3.5,   1.0,
  # Right Sock
  "Sock_BR",      5.5,   0.0,
  "Sock_BR",      6.5,   0.0,
  "Sock_BR",      6.5,   1.0,
  "Sock_BR",      5.5,   1.0
)

# Map region names to SampleLocation for merging with data
# (Left and right sides get the same value; _B suffix for back view)
region_to_location <- tribble(
  ~region,           ~SampleLocation, ~view,
  "Neck",            "Neck",          "Front",
  "Chest",           "Chest",         "Front",
  "Lower Chest",     "Lower Chest",   "Front",
  "Fly",             "Fly",           "Front",
  "Sleeve",          "Sleeve",        "Front",
  "Sleeve_R",        "Sleeve",        "Front",
  "Pant",            "Pant",          "Front",
  "Pant_R",          "Pant",          "Front",
  "Lower Pant",      "Lower Pant",    "Front",
  "Lower Pant_R",    "Lower Pant",    "Front",
  "Sock",            "Sock",          "Front",
  "Sock_R",          "Sock",          "Front",
  "Palm",            "Palm",          "Front",
  "Palm_R",          "Palm",          "Front",
  "Finger",          "Finger",        "Front",
  "Finger_R",        "Finger",        "Front",
  "Thumb",           "Thumb",         "Front",
  "Thumb_R",         "Thumb",         "Front",
  "Neck_B",          "Neck",          "Back",
  "Back",            "Back",          "Back",
  "Sleeve_B",        "Sleeve",        "Back",
  "Sleeve_BR",       "Sleeve",        "Back",
  "Pant_B",          "Pant",          "Back",
  "Pant_BR",         "Pant",          "Back",
  "Lower Pant_B",    "Lower Pant",    "Back",
  "Lower Pant_BR",   "Lower Pant",    "Back",
  "Sock_B",          "Sock",          "Back",
  "Sock_BR",         "Sock",          "Back"
)


# -----------------------------------------
# --- Step 2: Build plotting function    ---
# -----------------------------------------

# Helper: build a single paper doll plot (front or back)
plot_doll <- function(data_row, view = "Front", poly_df, region_map,
                      scale_limits, title_text = "") {

  # Filter geometry and mapping to this view
  if (view == "Front") {
    polys <- poly_df[["front"]]
  } else {
    polys <- poly_df[["back"]]
  }

  map_view <- region_map |> filter(view == !!view)

  # Join data values to regions
  polys_data <- polys |>
    left_join(map_view, by = "region") |>
    left_join(data_row, by = "SampleLocation")

  ggplot(polys_data, aes(x = x, y = y, group = region, fill = totalPAH_imputed)) +
    geom_polygon(color = "grey30", linewidth = 0.4) +
    scale_fill_gradientn(
      colors = c("white", "#FFFFCC", "#FFEDA0", "#FEB24C", "#FC4E2A", "#BD0026"),
      limits = scale_limits,
      na.value = "grey85",
      name = "Total PAH\n(Top 6, Imputed)"
    ) +
    labs(title = title_text) +
    coord_fixed() +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    )
}

# Package front and back geometry into a list
poly_list <- list(front = front_regions, back = back_regions)

# Helper: build a 2x2 page for one participant + one condition
# Layout: Front/Donning | Front/Doffing
#         Back/Donning  | Back/Doffing
plot_participant_page <- function(d_sub, participant, condition,
                                 poly_list, region_map, scale_limits) {

  d_pid <- d_sub |>
    filter(ParticipantID == participant, Condition == condition) |>
    select(SampleLocation, Timing_new, totalPAH_imputed)

  d_don <- d_pid |> filter(Timing_new == "Donning/Firefighting")
  d_dof <- d_pid |> filter(Timing_new == "Doffing")

  p_front_don <- plot_doll(d_don, "Front", poly_list, region_map, scale_limits,
                           "Front — Donning/FF")
  p_front_dof <- plot_doll(d_dof, "Front", poly_list, region_map, scale_limits,
                           "Front — Doffing")
  p_back_don  <- plot_doll(d_don, "Back", poly_list, region_map, scale_limits,
                           "Back — Donning/FF")
  p_back_dof  <- plot_doll(d_dof, "Back", poly_list, region_map, scale_limits,
                           "Back — Doffing")

  # Shared legend from one of the plots
  p_legend <- plot_doll(d_dof, "Front", poly_list, region_map, scale_limits, "") +
    theme(legend.position = "right")
  legend <- cowplot::get_legend(p_legend)

  combined <- (p_front_don + p_front_dof) / (p_back_don + p_back_dof) +
    plot_annotation(
      title    = paste0("Participant: ", participant),
      subtitle = paste0("Condition: ", condition,
                        "  |  Top 6 PAH (MLE-imputed)"),
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )

  combined
}

# -----------------------------------------
# --- Step 2b: Test with one participant ---
# -----------------------------------------

# Pick a participant with known detections for a good visual test
# BA24 had the highest exposure observations
test_plot <- plot_participant_page(
  d_sub      = d,
  participant = "BA24",
  condition   = "SS",
  poly_list   = poly_list,
  region_map  = region_to_location,
  scale_limits = c(global_min, global_max)
)

test_plot


# -----------------------------------------
# --- Step 2c: Updated plotting function ---
# -----------------------------------------

# Helper: compute centroids for text labels
get_centroids <- function(poly_df) {
  poly_df |>
    summarise(
      cx = mean(x),
      cy = mean(y),
      .by = region
    )
}

# Updated: build a single paper doll plot (front or back)
plot_doll <- function(data_row, view = "Front", poly_df, region_map,
                      scale_limits, title_text = "", show_legend = FALSE) {

  if (view == "Front") {
    polys <- poly_df[["front"]]
  } else {
    polys <- poly_df[["back"]]
  }

  map_view <- region_map |> filter(view == !!view)

  # Join: region → SampleLocation, then SampleLocation → data
  # Use many-to-many since left/right regions map to same SampleLocation
  polys_data <- polys |>
    left_join(map_view, by = "region", relationship = "many-to-many") |>
    left_join(data_row, by = "SampleLocation", relationship = "many-to-many")

  # Centroids for text labels
  centroids <- get_centroids(polys) |>
    left_join(map_view, by = "region", relationship = "many-to-many") |>
    left_join(data_row, by = "SampleLocation", relationship = "many-to-many") |>
    mutate(label = if_else(!is.na(totalPAH_imputed),
                           sprintf("%.2f", totalPAH_imputed), ""))

  p <- ggplot() +
    geom_polygon(data = polys_data,
                 aes(x = x, y = y, group = region, fill = totalPAH_imputed),
                 color = "grey30", linewidth = 0.4) +
    geom_text(data = centroids,
              aes(x = cx, y = cy, label = label),
              size = 2.2, color = "black", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("white", "#FFFFCC", "#FFEDA0", "#FEB24C", "#FC4E2A", "#BD0026"),
      limits = scale_limits,
      na.value = "grey85",
      name = "Total PAH\n(Top 6, Imputed)"
    ) +
    labs(title = title_text) +
    coord_fixed() +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}

# Updated: build a 2x2 page for one participant + one condition
plot_participant_page <- function(d_sub, participant, condition,
                                 poly_list, region_map, scale_limits) {

  d_pid <- d_sub |>
    filter(ParticipantID == participant, Condition == condition) |>
    select(SampleLocation, Timing_new, totalPAH_imputed)

  d_don <- d_pid |> filter(Timing_new == "Donning/Firefighting")
  d_dof <- d_pid |> filter(Timing_new == "Doffing")

  p_front_don <- plot_doll(d_don, "Front", poly_list, region_map, scale_limits,
                           "Front — Donning/FF")
  p_front_dof <- plot_doll(d_dof, "Front", poly_list, region_map, scale_limits,
                           "Front — Doffing")
  p_back_don  <- plot_doll(d_don, "Back", poly_list, region_map, scale_limits,
                           "Back — Donning/FF")
  p_back_dof  <- plot_doll(d_dof, "Back", poly_list, region_map, scale_limits,
                           "Back — Doffing", show_legend = TRUE)

  combined <- (p_front_don + p_front_dof) / (p_back_don + p_back_dof) +
    plot_annotation(
      title    = paste0("Participant: ", participant),
      subtitle = paste0("Condition: ", condition,
                        "  |  Top 6 PAH (MLE-imputed)"),
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )

  combined
}

# Re-test with BA24/SS
test_plot <- plot_participant_page(
  d_sub        = d,
  participant  = "BA24",
  condition    = "SS",
  poly_list    = poly_list,
  region_map   = region_to_location,
  scale_limits = c(global_min, global_max)
)

test_plot

d |>
  filter(ParticipantID == "BA24", Condition == "SS") |>
  select(SampleLocation, Timing_new, totalPAH_imputed) |>
  arrange(SampleLocation, Timing_new)


# -----------------------------------------
# --- Step 2d: Fix back view mapping    ---
# -----------------------------------------

# Update region_to_location: only Back maps to data on the back view
# All other back regions remain in geometry but get no SampleLocation → fill = NA → grey
region_to_location <- tribble(
  ~region,           ~SampleLocation, ~view,
  # Front view mappings
  "Neck",            "Neck",          "Front",
  "Chest",           "Chest",         "Front",
  "Lower Chest",     "Lower Chest",   "Front",
  "Fly",             "Fly",           "Front",
  "Sleeve",          "Sleeve",        "Front",
  "Sleeve_R",        "Sleeve",        "Front",
  "Pant",            "Pant",          "Front",
  "Pant_R",          "Pant",          "Front",
  "Lower Pant",      "Lower Pant",    "Front",
  "Lower Pant_R",    "Lower Pant",    "Front",
  "Sock",            "Sock",          "Front",
  "Sock_R",          "Sock",          "Front",
  "Palm",            "Palm",          "Front",
  "Palm_R",          "Palm",          "Front",
  "Finger",          "Finger",        "Front",
  "Finger_R",        "Finger",        "Front",
  "Thumb",           "Thumb",         "Front",
  "Thumb_R",         "Thumb",         "Front",
  # Back view — only Back maps to data
  "Back",            "Back",          "Back"
)

# Re-test
test_plot <- plot_participant_page(
  d_sub        = d,
  participant  = "BA24",
  condition    = "SS",
  poly_list    = poly_list,
  region_map   = region_to_location,
  scale_limits = c(global_min, global_max)
)

test_plot
