
# Prepare data - combine Shirt with Sleeve, use new_totalPAH
pah_by_location <- qc %>%
  mutate(`Sample Location` = ifelse(`Sample Location` == "Shirt", "Sleeve", `Sample Location`)) %>%
  group_by(`Sample Location`) %>%
  summarise(new_totalPAH = sum(new_totalPAH, na.rm = TRUE), .groups = "drop")

# ===== FRONT VIEW REGIONS =====
body_regions_front <- tribble(
  ~region_name, ~x, ~y, ~width, ~height,
  "Neck", 0, 9, 0.8, 0.6,
  "Chest", 0, 6.75, 2.5, 3,        
  "Lower Chest", 0, 4.5, 2.5, 1.5, 
  "Sleeve", -2.5, 5, 1.2, 3,
  "Palm", -2.5, 2.5, 0.8, 1,
  "Finger", -2.5, 1.5, 0.8, 0.8,
  "Thumb", -3.2, 2.5, 0.5, 0.6,
  "Pant", -0.8, 2, 1.2, 2,
  "Lower Pant", -0.8, 0, 1.2, 1.5,
  "Sock", -0.8, -1, 1.2, 0.8,
  "Sleeve", 2.5, 5, 1.2, 3,
  "Palm", 2.5, 2.5, 0.8, 1,
  "Finger", 2.5, 1.5, 0.8, 0.8,
  "Thumb", 3.2, 2.5, 0.5, 0.6,
  "Pant", 0.8, 2, 1.2, 2,
  "Lower Pant", 0.8, 0, 1.2, 1.5,
  "Sock", 0.8, -1, 1.2, 0.8,
  "Fly", 0, 3.25, 1.5, 0.5
)

# ===== BACK VIEW REGIONS =====
body_regions_back <- tribble(
  ~region_name, ~x, ~y, ~width, ~height,
  "Neck", 0, 9, 0.8, 0.6,
  "Back", 0, 6, 2.5, 4.5,
  "Sleeve", -2.5, 5, 1.2, 3,
  "Palm", -2.5, 2.5, 0.8, 1,
  "Finger", -2.5, 1.5, 0.8, 0.8,
  "Thumb", -3.2, 2.5, 0.5, 0.6,
  "Pant", -0.8, 2, 1.2, 2,
  "Lower Pant", -0.8, 0, 1.2, 1.5,
  "Sock", -0.8, -1, 1.2, 0.8,
  "Sleeve", 2.5, 5, 1.2, 3,
  "Palm", 2.5, 2.5, 0.8, 1,
  "Finger", 2.5, 1.5, 0.8, 0.8,
  "Thumb", 3.2, 2.5, 0.5, 0.6,
  "Pant", 0.8, 2, 1.2, 2,
  "Lower Pant", 0.8, 0, 1.2, 1.5,
  "Sock", 0.8, -1, 1.2, 0.8,
  "Fly", 0, 3.25, 1.5, 0.5
)

# ===== FRONT PLOT =====
plot_data_front <- body_regions_front %>%
  left_join(pah_by_location, by = c("region_name" = "Sample Location")) %>%
  replace_na(list(new_totalPAH = 0)) %>%
  group_by(x, y, width, height) %>%
  summarise(
    region_name = first(region_name),
    new_totalPAH = max(new_totalPAH),
    .groups = "drop"
  )

# ===== BACK PLOT =====
plot_data_back <- body_regions_back %>%
  left_join(pah_by_location, by = c("region_name" = "Sample Location")) %>%
  mutate(
    new_totalPAH = ifelse(region_name == "Back", new_totalPAH, NA),
    label_text = ifelse(region_name == "Back", 
                        paste0(region_name, "\n", round(new_totalPAH, 2)), 
                        "")
  ) %>%
  replace_na(list(new_totalPAH = 0))

# Calculate shared color scale limits
max_pah <- max(plot_data_front$new_totalPAH, na.rm = TRUE)

# Create shared fill scale
shared_scale <- scale_fill_gradientn(
  colors = c("white", "yellow", "orange", "red", "darkred"),
  values = scales::rescale(c(0, 0.5, 2, 5, max_pah)),
  limits = c(0, max_pah),
  name = "Total PAH\n(µg/g)",
  na.value = "white"
)

front_plot <- ggplot(plot_data_front, aes(x = x, y = y, fill = new_totalPAH)) +
  geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
  geom_text(aes(label = paste0(region_name, "\n", round(new_totalPAH, 2))), 
            size = 2.5, lineheight = 0.8) +
  shared_scale +
  coord_fixed() +
  theme_void() +
  labs(title = "Front View: Total PAH Exposure")

back_plot <- ggplot(plot_data_back, aes(x = x, y = y, fill = new_totalPAH)) +
  geom_tile(aes(width = width, height = height), color = "black", linewidth = 0.5) +
  geom_text(aes(label = label_text), size = 2.5, lineheight = 0.8) +
  shared_scale +
  guides(fill = "none") +  # Remove legend from back plot
  coord_fixed() +
  theme_void() +
  labs(title = "Back View: Total PAH Exposure")

# Display both plots - legend will be between them

front_plot + back_plot