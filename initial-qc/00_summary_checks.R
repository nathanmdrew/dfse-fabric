# install.packages("tidyverse")
# install.packages("gt")
# install.packages("writexl")
#install.packages("emmeans")
#install.packages("car")
#install.packages("rstatix")
#install.packages("ggpubr")

library(dplyr)
library(readxl)
library(purrr)
library(gt)
library(writexl)
library(broom)
library(emmeans)
library(car)
library(rstatix)
library(ggpubr)
library(tidyr)
library(stringr)

data_fabric_raw <- read_excel(path="S:/NIOSH/DFSE Project/Human_Subject_Results_Ver5_master.xlsx",
                       sheet="Fabric",
                       col_names = TRUE)

temp_summary <- data_fabric_raw %>% group_by(SampleID) %>% summarize(num=n())

pah1_summary <- data_fabric_raw %>% summarize(n=n(),
                                              mean=mean(Acenaphthene...23),
                                              `std dev`=sqrt(var(Acenaphthene...23)),
                                              median=median(Acenaphthene...23),
                                              minimum=min(Acenaphthene...23),
                                              max=max(Acenaphthene...23)
                                              )
str(data_fabric_raw)

summary_table <- map_dfr(names(data_fabric_raw)[23:37], function(col) {
  x <- data_fabric_raw[[col]]
  tibble(
    Analyte = col,
    n       = sum(!is.na(x)),
    n_missing  = sum(is.na(x)),
    mean    = mean(x, na.rm = TRUE),
    sd      = sd(x, na.rm = TRUE),
    median  = median(x, na.rm = TRUE),
    min     = min(x, na.rm = TRUE),
    max     = max(x, na.rm = TRUE)
  )
})

summary_table

gt_table <- summary_table %>%
  gt() %>%
  tab_header(
    title = "Summary Statistics of Analytes",
    subtitle = "Exploratory Data Analysis (Columns 23–37), units in μg/g"
  ) %>%
  fmt_number(
    columns = c(mean, sd, median, min, max),
    decimals = 4
  ) %>%
  cols_label(
    Analyte   = "Analyte",
    n         = "Non-Missing",
    n_missing = "Missing",
    mean      = "Mean",
    sd        = "Std Dev",
    median    = "Median",
    min       = "Minimum",
    max       = "Maximum"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

gt_table

# Export the raw summary_table data to Excel
write_xlsx(summary_table, "S:\\NIOSH\\DFSE Project\\summary_table.xlsx")
# Save to PDF
gtsave(gt_table, "S:\\NIOSH\\DFSE Project\\summary_table.pdf")


# Factor-ize the group variables
data_fabric_raw$`Sample Location` <- as.factor(data_fabric_raw$`Sample Location`)
data_fabric_raw$Condition <- as.factor(data_fabric_raw$Condition)
data_fabric_raw$Timing <- as.factor(data_fabric_raw$Timing)

summary_table2 <- data_fabric_raw %>%
  group_by(Condition, `Sample Location`, Timing) %>%
  summarise(
    N           = n(),
    N_non_detect= sum(`Total_PAH_ug/g` == 0, na.rm = TRUE),
    mean        = mean(`Total_PAH_ug/g`, na.rm = TRUE),
    sd          = sd(`Total_PAH_ug/g`, na.rm = TRUE),
    median      = median(`Total_PAH_ug/g`, na.rm = TRUE),
    min         = min(`Total_PAH_ug/g`, na.rm = TRUE),
    max         = max(`Total_PAH_ug/g`, na.rm = TRUE),
    .groups = "drop"
  )

summary_table2

gt_table2 <- summary_table2 %>%
  gt() %>%
  tab_header(
    title = "Summary of Total PAH (μg/g) by Condition, Location, and Timing"
  ) %>%
  fmt_number(
    columns = c(mean, sd, median, min, max),
    decimals = 2
  ) %>%
  cols_label(
    Condition         = "Condition",
    `Sample Location` = "Sample Location",
    Timing            = "Timing",
    N                 = "N",
    N_non_detect      = "N Non-Detect",
    mean              = "Mean",
    sd                = "Std Dev",
    median            = "Median",
    min               = "Minimum",
    max               = "Maximum"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

gt_table2

# Export the raw summary_table data to Excel
write_xlsx(summary_table2, "S:\\NIOSH\\DFSE Project\\summary_table2_totalPAH_condition_location_timing.xlsx")
# Save to PDF
gtsave(gt_table2, "S:\\NIOSH\\DFSE Project\\summary_table2.pdf")




#### research question 1
# Fit linear model
anova_model <- lm(`Total_PAH_ug/g` ~ Condition, data = data_fabric_raw)

# Run ANOVA
anova_results <- anova(anova_model)
anova_results #p-value 0.1463, no sig differences

# Tukey HSD test for pairwise differences
tukey_results <- TukeyHSD(aov(anova_model))
tukey_results

# ANOVA results (tidy format)
tidy(anova(anova_model))

# Pairwise comparisons with Tukey adjustment
emmeans(anova_model, pairwise ~ Condition, adjust = "tukey")

# Normality check (QQ plot)
plot(anova_model, which = 2)

# Equal variance test
leveneTest(`Total_PAH_ug/g` ~ Condition, data = data_fabric_raw)

# Kruskal-Wallis test
kruskal_results <- data_fabric_raw %>%
  kruskal_test(`Total_PAH_ug/g` ~ Condition)
kruskal_results

# Pairwise Wilcoxon rank-sum tests with Benjamini-Hochberg adjustment
pairwise_results <- data_fabric_raw %>%
  pairwise_wilcox_test(`Total_PAH_ug/g` ~ Condition, p.adjust.method = "BH")
pairwise_results

# Format results for plotting
pairwise_results <- pairwise_results %>%
  add_xy_position(x = "Condition")

# Plot boxplot with stats
ggplot(data_fabric_raw, aes(x = Condition, y = `Total_PAH_ug/g`, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black") +
  stat_pvalue_manual(pairwise_results, hide.ns = TRUE, tip.length = 0.01) +
  labs(
    title = "Total PAH (μg/g) by Condition",
    subtitle = paste0("Kruskal-Wallis p = ", signif(kruskal_results$p, 3)),
    x = "Condition",
    y = "Total PAH"
  ) +
  theme_minimal() +
  theme(legend.position = "none")





#### research question 2
# Fit linear model
anova_model <- lm(`Total_PAH_ug/g` ~ Sample_Location, data = data_fabric_raw)

# Run ANOVA
anova_results <- anova(anova_model) #p-value 0.000921

# ANOVA results (tidy format)
tidy(anova(anova_model))

# Pairwise comparisons with Tukey adjustment
emmeans(anova_model, pairwise ~ `Sample Location`, adjust = "tukey")

# Normality check (QQ plot)
plot(anova_model, which = 2)

# Kruskal-Wallis test
kruskal_results <- data_fabric_raw %>%
  kruskal_test(`Total_PAH_ug/g` ~ Sample_Location)
kruskal_results

### see program 00b




#### research question 3
subset <- data_fabric_raw %>% filter(Sample_Location %in% c("Finger", "Palm", "Thumb"))

# Fit linear model
anova_model <- lm(`Total_PAH_ug/g` ~ Timing, data = data_fabric_raw)

# ANOVA results (tidy format)
tidy(anova(anova_model))

### see program 00c








# make column names safe (replace spaces with underscores)
names(data_fabric_raw) <- gsub(" ", "_", names(data_fabric_raw))

# select the analyte columns (23:37 = your 15 numeric analytes)
analytes <- names(data_fabric_raw)[23:37]

summary_table <- data_fabric_raw %>%
  select(Sample_Location, all_of(analytes)) %>%
  pivot_longer(cols = all_of(analytes), names_to = "Analyte", values_to = "Value") %>%
  group_by(Sample_Location, Analyte) %>%
  summarise(
    N      = sum(!is.na(Value)),
    mean   = mean(Value, na.rm = TRUE),
    sd     = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    min    = min(Value, na.rm = TRUE),
    max    = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Sample_Location, Analyte)

# Export the raw summary_table data to Excel
write_xlsx(summary_table, "S:\\NIOSH\\DFSE Project\\summary_table_pah_by_sample_location.xlsx")
                                     