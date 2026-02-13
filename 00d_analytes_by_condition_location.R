data_fabric_raw %>% distinct(Sample_Location)

# Replace spaces in column names if needed
names(data_fabric_raw) <- gsub(" ", "_", names(data_fabric_raw))

# Suppose your analytes are columns 23:37
analytes <- names(data_fabric_raw)[23:37]

summary_stats <- data_fabric_raw %>%
  select(Condition, Sample_Location, all_of(analytes)) %>%
  pivot_longer(cols = all_of(analytes), names_to = "Analyte", values_to = "Value") %>%
  group_by(Sample_Location, Condition, Analyte) %>%
  summarise(
    N = sum(!is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Sample_Location, desc(mean))

top8_analytes <- summary_stats %>%
  group_by(Sample_Location, Analyte) %>%
  summarise(mean_overall = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample_Location) %>%
  slice_max(mean_overall, n = 8) %>%
  pull(Analyte)

summary_top8 <- summary_stats %>%
  filter(Analyte %in% top8_analytes)

top8_per_combo <- summary_stats %>%
  group_by(Sample_Location, Condition) %>%
  slice_max(order_by = mean, n = 8) %>%   # top 8 by mean
  summarise(
    top_analytes = str_c(Analyte, collapse = "\n"),  # newline separated
    .groups = "drop"
  )

ggplot(top8_per_combo, aes(x = Condition, y = Sample_Location, label = top_analytes)) +
  geom_text(size = 3, hjust = 0, vjust = 0.5, lineheight = 0.9) +
  scale_y_discrete(limits = rev(unique(top8_per_combo$Sample_Location))) + # optional: reverse row order
  labs(
    x = "Condition",
    y = "Sample Location",
    title = "Top 8 Analytes per Sample Location × Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )







top8_per_combo <- summary_stats %>%
  group_by(Sample_Location, Condition) %>%
  slice_max(order_by = mean, n = 8) %>%
  arrange(desc(mean)) %>%  # optional: sort by mean
  summarise(
    top_analytes = str_c(Analyte, collapse = "\n"),  # newline-separated
    .groups = "drop"
  )

ggplot(top8_per_combo, aes(x = Condition, y = Sample_Location)) +
  geom_tile(fill = "#f0f0f0", color = "white", width = 0.9, height = 0.9) +
  geom_text(aes(label = top_analytes), size = 3, lineheight = 0.9, hjust = 0) +
  scale_y_discrete(limits = rev(unique(top8_per_combo$Sample_Location))) + # reverse rows
  labs(
    x = "Condition",
    y = "Sample Location",
    title = "Top 8 Analytes per Sample Location × Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    plot.title = element_text(face = "bold", size = 14)
  )



pdf("S:\\NIOSH\\DFSE Project\\Top8_Analytes_by_SampleLocation_Condition.pdf",
    width = 12, height = 30)

# Print the plot (same as above)
print(
  ggplot(top8_per_combo, aes(x = Condition, y = Sample_Location)) +
    geom_tile(fill = "#f0f0f0", color = "white", width = 0.9, height = 0.9) +
    geom_text(aes(label = top_analytes), size = 3, lineheight = 0.9, hjust = 0) +
    scale_y_discrete(limits = rev(unique(top8_per_combo$Sample_Location))) + # reverse rows
    labs(
      x = "Condition",
      y = "Sample Location",
      title = "Top 8 Analytes per Sample Location × Condition"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA),
      plot.title = element_text(face = "bold", size = 14)
    )
)

dev.off()




# strictly 8
top8_per_combo <- summary_stats %>%
  group_by(Sample_Location, Condition) %>%
  arrange(desc(mean), Analyte) %>%  # sort by mean descending, then analyte name
  slice_head(n = 8) %>%             # exactly 8 rows per group
  summarise(
    top_analytes = str_c(Analyte, collapse = "\n"),
    top_mean = mean(mean, na.rm = TRUE),
    .groups = "drop"
  )

pdf("S:\\NIOSH\\DFSE Project\\Top8_Analytes_by_SampleLocation_Condition_STRICT.pdf",
    width = 12, height = 30)



print(
  ggplot(top8_per_combo, aes(x = Condition, y = Sample_Location, fill = top_mean)) +
    geom_tile(color = "white", width = 0.9, height = 0.9) +
    geom_text(aes(label = top_analytes), size = 3, lineheight = 0.9, hjust = 0) +
    scale_fill_viridis_c(option = "C") +
    scale_y_discrete(limits = rev(unique(top8_per_combo$Sample_Location))) +
    labs(
      x = "Condition",
      y = "Sample Location",
      fill = "Mean of Top 8",
      title = "Top 8 Analytes per Sample Location × Condition"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )
)

dev.off()
