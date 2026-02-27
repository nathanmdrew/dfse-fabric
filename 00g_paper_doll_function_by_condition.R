library(tidyverse)
library(patchwork)

# ============================================================
# HELPER: Geometric mean (excluding zeros/non-detects)
# Assumption: zeros represent non-detects and are excluded
# from the geometric mean calculation. If all values for a
# location/condition are zero, geomean returns 0.
# ============================================================
geom_mean <- function(x, na.rm = TRUE) {
  x_pos <- x[x > 0 & !is.na(x)]
  if (length(x_pos) == 0) return(0)
  exp(mean(log(x_pos), na.rm = na.rm))
}

# ============================================================
# SUMMARIZE: Geometric mean by Condition and Sample Location
# ============================================================
pah_condition_summary <- d %>%
  mutate(`Sample Location` = ifelse(`Sample Location` == "Shirt", "Sleeve", `Sample Location`)) %>%
  group_by(Condition, `Sample Location`) %>%
  summarise(
    geomean_PAH = geom_mean(new_totalPAH),
    n_subjects  = n_distinct(ParticipantID),
    n_nonzero   = sum(new_totalPAH > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Quick check
print(pah_condition_summary)

# ============================================================
# FUNCTION: build_condition_paper_doll
# Takes a single-condition summary and returns a patchwork plot
# ============================================================
build_condition_paper_doll <- function(condition_summary, condition) {
  
  # ===== FRONT VIEW REGIONS =====
  body_regions_front <- tribble(
    ~region_name, ~x,   ~y,    ~width, ~height,
    "Neck",        0,    9,     0.8,    0.6,
    "Chest",       0,    6.75,  2.5,    3,
    "Lower Chest", 0,    4.5,   2.5,    1.5,
    "Sleeve",     -2.5,  5,     1.2,    3,
    "Palm",       -2.5,  2.5,   0.8,    1,
    "Finger",     -2.5,  1.5,   0.8,    0.8,
    "Thumb",      -3.2,  2.5,   0.5,    0.6,
    "Pant",       -0.8,  2,     1.2,    2,
    "Lower Pant", -0.8,  0,     1.2,    1.5,
    "Sock",       -0.8, -1,     1.2,    0.8,
    "Sleeve",      2.5,  5,     1.2,    3,
    "Palm",        2.5,  2.5,   0.8,    1,
    "Finger",      2.5,  1.5,   0.8,    0.8,
    "Thumb",       3.2,  2.5,   0.5,    0.6,
    "Pant",        0.8,  2,     1.2,    2,
    "Lower Pant",  0.8,  0,     1.2,    1.5,
    "Sock",        0.8, -1,     1.2,    0.8,
    "Fly",         0,    3.25,  1.5,    0.5
  )
  
  # ===== BACK VIEW REGIONS =====
  body_regions_back <- tribble(
    ~region_name, ~x,   ~y,    ~width, ~height,
    "Neck",        0,    9,     0.8,    0.6,
    "Back",        0,    6,     2.5,    4.5,
    "Sleeve",     -2.5,  5,     1.2,    3,
    "Palm",       -2.5,  2.5,   0.8,    1,
    "Finger",     -2.5,  1.5,   0.8,    0.8,
    "Thumb",      -3.2,  2.5,   0.5,    0.6,
    "Pant",       -0.8,  2,     1.2,    2,
    "Lower Pant", -0.8,  0,     1.2,    1.5,
    "Sock",       -0.8, -1,     1.2,    0.8,
    "Sleeve",      2.5,  5,     1.2,    3,
    "Palm",        2.5,  2.5,   0.8,    1,
    "Finger",      2.5,  1.5,   0.8,    0.8,
    "Thumb",       3.2,  2.5,   0.5,    0.6,
    "Pant",        0.8,  2,     1.2,    2,
    "Lower Pant",  0.8,  0,     1.2,    1.5,
    "Sock",        0.8, -1,     1.2,    0.8,
    "Fly",         0,    3.25,  1.5,    0.5
  )
  
  # ===== PREPARE PLOT DATA =====
  plot_data_front <- body_regions_front %>%
    left_join(condition_summary, by = c("region_name" = "Sample Location")) %>%
    replace_na(list(geomean_PAH = 0)) %>%
    group_by(x, y, width, height) %>%
    summarise(
      region_name = first(region_name),
      geomean_PAH = max(geomean_PAH),
      .groups = "drop"
    )
  
  plot_data_back <- body_regions_back %>%
    left_join(condition_summary, by = c("region_name" = "Sample Location")) %>%
    mutate(
      geomean_PAH = ifelse(region_name == "Back", geomean_PAH, NA),
      label_text  = ifelse(region_name == "Back",
                           paste0(region_name, "\n", round(geomean_PAH, 2)),
                           "")
    ) %>%
    replace_na(list(geomean_PAH = 0))
  
  return(list(front = plot_data_front, back = plot_data_back))
}

# ============================================================
# BUILD ALL 3 CONDITIONS AND USE A SHARED COLOR SCALE
# ============================================================
conditions <- c("OL", "SS", "SL")

# Build plot data for all conditions first so we can get a
# global max for a consistent color scale across all 3 plots
all_plot_data <- map(conditions, function(cond) {
  cond_summary <- pah_condition_summary %>% filter(Condition == cond)
  build_condition_paper_doll(cond_summary, cond)
}) %>%
  set_names(conditions)

# Calculate global max across all conditions for shared scale
global_max <- map_dbl(all_plot_data, ~ max(.x$front$geomean_PAH, na.rm = TRUE)) %>% max()
global_max <- max(global_max, 1)  # Floor of 1 to avoid degenerate scale

# Shared color scale
shared_scale <- scale_fill_gradientn(
  colors   = c("white", "yellow", "orange", "red", "darkred"),
  values   = scales::rescale(c(0, 0.5, 2, 5, global_max)),
  limits   = c(0, global_max),
  name     = "Geom. Mean\nTotal PAH\n(µg/g)",
  na.value = "white"
)

# ============================================================
# BUILD PLOTS FOR EACH CONDITION
# ============================================================
condition_plots <- map(conditions, function(cond) {
  
  pd_front <- all_plot_data[[cond]]$front
  pd_back  <- all_plot_data[[cond]]$back
  
  front_plot <- ggplot(pd_front, aes(x = x, y = y, fill = geomean_PAH)) +
    geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
    geom_text(aes(label = paste0(region_name, "\n", round(geomean_PAH, 2))),
              size = 2.5, lineheight = 0.8) +
    shared_scale +
    coord_fixed() +
    theme_void() +
    labs(subtitle = "Front View")
  
  back_plot <- ggplot(pd_back, aes(x = x, y = y, fill = geomean_PAH)) +
    geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
    geom_text(aes(label = label_text), size = 2.5, lineheight = 0.8) +
    shared_scale +
    guides(fill = "none") +
    coord_fixed() +
    theme_void() +
    labs(subtitle = "Back View")
  
  front_plot + back_plot +
    plot_annotation(
      title = paste0("Typical Subject | Condition: ", cond,
                     "  (Geometric Mean Total PAH, n=", 
                     n_distinct(d$ParticipantID[d$Condition == cond]), " subjects)"),
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
    )
})

# ============================================================
# VIEW INDIVIDUAL PLOTS
# ============================================================
condition_plots[["OL"]]
condition_plots[["SS"]]
condition_plots[["SL"]]

# ============================================================
# SAVE TO A SINGLE PDF (one page per condition)
# ============================================================
pdf("PAH_paper_doll_by_condition.pdf", width = 10, height = 8)
walk(condition_plots, print)
dev.off()

message("Saved: PAH_paper_doll_by_condition.pdf")