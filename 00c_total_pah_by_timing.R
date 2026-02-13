

# --- 1. Identify valid groups ---
group_summary <- data_fabric_raw %>%
  group_by(Timing) %>%
  summarise(
    n_non_missing = sum(!is.na(`Total_PAH_ug/g`)),
    n_unique = n_distinct(`Total_PAH_ug/g`, na.rm = TRUE),
    .groups = "drop"
  )

valid_conditions <- group_summary %>%
  filter(n_non_missing >= 2, n_unique >= 2) %>%
  pull(Timing)

# --- 2. Filter dataset ---
filtered_data <- data_fabric_raw %>%
  filter(Timing %in% valid_conditions,
         Timing != "NA")

# --- 3. Run Kruskal-Wallis ---
kruskal_res <- filtered_data %>%
  kruskal_test(`Total_PAH_ug/g` ~ Timing)

# --- 4. Prepare all pairwise comparisons ---
condition_pairs <- combn(valid_conditions, 2, simplify = FALSE)

# Check each pair for enough data & variability
valid_pairs <- map_dfr(condition_pairs, function(pair) {
  x <- na.omit(filtered_data$`Total_PAH_ug/g`[filtered_data$Timing == pair[1]])
  y <- na.omit(filtered_data$`Total_PAH_ug/g`[filtered_data$Timing == pair[2]])
  
  if(length(x) >= 2 & length(y) >= 2 & length(unique(x)) > 1 & length(unique(y)) > 1) {
    tibble(group1 = pair[1], group2 = pair[2])
  } else {
    NULL  # exclude invalid pairs
  }
})

# --- 5. Run Wilcoxon for valid pairs ---
pairwise_res <- map_dfr(1:nrow(valid_pairs), function(i) {
  x <- na.omit(filtered_data$`Total_PAH_ug/g`[filtered_data$Timing == valid_pairs$group1[i]])
  y <- na.omit(filtered_data$`Total_PAH_ug/g`[filtered_data$Timing == valid_pairs$group2[i]])
  res <- broom::tidy(wilcox.test(x, y))
  tibble(
    group1 = valid_pairs$group1[i],
    group2 = valid_pairs$group2[i],
    statistic = res$statistic,
    p.value = res$p.value
  )
}) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

# --- 6. Prepare y-position for plotting ---
y_max <- filtered_data %>%
  group_by(Timing) %>%
  summarise(ymax = max(`Total_PAH_ug/g`, na.rm = TRUE), .groups = "drop")

plot_pairwise <- pairwise_res %>%
  left_join(y_max, by = c("group1" = "Timing")) %>%
  rename(y1 = ymax) %>%
  left_join(y_max, by = c("group2" = "Timing")) %>%
  rename(y2 = ymax) %>%
  mutate(y.position = pmax(y1, y2) * 1.05) %>%
  select(group1, group2, p.adj, y.position)

# Round p-values to 3 decimal places
plot_pairwise <- plot_pairwise %>%
  mutate(
    p.label = ifelse(p.adj < 0.001, "<0.001", sprintf("%.3f", p.adj))
  )

# --- 7. Plot boxplot with means and significance bars ---
ggplot(filtered_data, aes(x = Timing, y = `Total_PAH_ug/g`, fill = Timing)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black") +
  ggpubr::stat_pvalue_manual(plot_pairwise, label = "p.label", hide.ns = TRUE,
                             tip.length=0.01, step.increase=0.05) +
  labs(
    title = "Total PAH by Timing",
    subtitle = paste0("Kruskal-Wallis p = ", signif(kruskal_res$p, 3)),
    x = "Timing",
    y = "Total PAH"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
