### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-02-13
### Purpose: Continue EDA of the firefighter fabric data.
###          Previous EDA treated each point as independent, and did not reflect
###          the grouped design (n=23 subjects). Information was received to
###          identify individual subjects, and some variable clean-up was
###          identified (viz. 00_eda.R).
###          
###          This analysis will revisit those completed previously while also
###          utilzing the subjects [viz. 00b, 00c, 00d R analyses]
###
### TODOs:   LoDs - what does "below LoD" mean? Per observation or a central
###          measure?
###        

library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
library(ggplot2)

setwd("C:/Users/vom8/dfse-fabric/")

### let's get a clean copy of the fabric data to always work from
#original_data <- read_excel(path="C:/Users/vom8/dfse-fabric/Human_Subject_Results_Ver5_master.xlsx",
#                            sheet="Fabric")

#saveRDS(original_data, file="original_data.RDS")
#rm(original_data)

# Parse SampleID into components (keep SampleID as the original full value)
d <- readRDS("original_data.RDS")
d$rownum <- seq(from=1, to=nrow(d))

d <- d %>%
  mutate(SampleID = as.character(SampleID)) %>%
  mutate(ParticipantID = str_extract(SampleID, "^[^_]+")) %>%
  mutate(rest = str_remove(SampleID, "^[^_]+_")) %>%
  separate(rest, into = c("Condition2", "SampleLocation2", "SampleID2_raw"),
           sep = "_", remove = TRUE, extra = "merge", fill = "right") %>%
  mutate(SampleID2 = as.integer(str_extract(SampleID2_raw, "\\d+"))) %>%
  select(-SampleID2_raw)

# Fix ParticipantID - uppercase, and fix the O instead of a 0.
d <- d %>%
  mutate(ParticipantID = toupper(ParticipantID)) %>%
  mutate(ParticipantID = str_replace_all(ParticipantID, "O", "0"))

# Already checked that SampleLocation2 and Condition2 match the existing variables
# but should allow for a separation for the "Shirt" occurrence of S
d <- d %>% select(-Condition2, -SampleLocation2)

# Recalculate total PAH; 5 rows where totals appear to be incorrect
d <- d %>% mutate(new_totalPAH = rowSums(select(., 23:37), na.rm = FALSE))

# drop the column with the single value of 1.1
d <- d %>% select(-`...54`)


################################
### QC Stuff
################################
#qc <- unique(d$ParticipantID)
# Good, 23 subjects

qc_table <- table(d$`Sample Location`)
qc_table <- arrange(qc_table)
# Note fly=8; lower chest=2; shirt=1

# Check if Total_PAH_ug/g equals the sum of columns 23-37
d_qc <- d %>%
  mutate(calculated_total = rowSums(select(., 23:37), na.rm = FALSE))

# View mismatches and NAs
qc_total <- d_qc %>%
  mutate(match = `Total_PAH_ug/g` == calculated_total) %>%
  filter(!match | is.na(`Total_PAH_ug/g`)) %>%
  select(rownum, ParticipantID, SampleID, `Total_PAH_ug/g`, calculated_total, match)

qc_total <- qc_total %>% filter(rownum %in% c(334,353,379,419,442))
qc_total$rownum <- qc_total$rownum+1
write.csv(qc_total, file="QC_totalPAH_20260217.csv")
# most actually match - check out rows 334, 353, 379, 419 (***), 442
# row 419 is the only one thats not currently missing, but the sum is different.
# the other 4 are missing the total.
# note the Excel rows are +1 due to the header row




# 1. Histogram of Total_PAH_ug/g across all samples
summary(d$`Total_PAH_ug/g`) #4 missing
qc <- d %>% filter(is.na(`Total_PAH_ug/g`)) # 4 missing
names(d)

ggplot(d, aes(x = `Total_PAH_ug/g`)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
  stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), 
           vjust = -0.5, size = 3.5) +
  labs(title = "Distribution of Total PAH",
       subtitle = "Bins of 0.5 µg/g",
       x = "Total PAH (µg/g)",
       y = "Frequency") +
  theme_minimal()

# Create histogram object to extract bin information
h <- hist(d$`Total_PAH_ug/g`, breaks = 20, plot = FALSE)

# Create a summary table
bin_summary <- tibble(
  bin_min = h$breaks[1:(length(h$breaks)-1)],
  bin_max = h$breaks[2:length(h$breaks)],
  count = h$counts,
  frequency = h$counts / sum(h$counts, na.rm = TRUE)
)

print(bin_summary)

# 2. Histograms for each ParticipantID (faceted)
ggplot(d, aes(x = `Total_PAH_ug/g`)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "black", alpha = 0.7) +
  stat_bin(bins = 10, geom = "text", aes(label = after_stat(count)), 
           vjust = -0.5, size = 2.5) +
  facet_wrap(~ParticipantID, ncol = 5) +
  labs(title = "Distribution of Total PAH by Subject",
       x = "Total PAH (µg/g)",
       y = "Frequency") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8))

# 3. Histogram of Total_PAH_ug/g across all samples on log scale (zeros excluded)
ggplot(d, aes(x = `Total_PAH_ug/g`)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
  stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), 
           vjust = -0.5, size = 3.5) +
  scale_x_log10() +
  labs(title = "Distribution of Total PAH (log scale)",
       x = "Total PAH (µg/g, log scale)",
       y = "Frequency") +
  theme_minimal()

# 4. Histogram of Total_PAH_ug/g across all samples on asinh scale
ggplot(d, aes(x = `Total_PAH_ug/g`)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
  stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), 
           vjust = -0.5, size = 3.5) +
  scale_x_continuous(trans = scales::asinh_trans()) +
  labs(title = "Distribution of Total PAH (asinh scale)",
       x = "Total PAH (µg/g, asinh scale)",
       y = "Frequency") +
  theme_minimal()


# 5. Histograms of each PAH

pah_cols <- names(d)[23:37]

for (pah in pah_cols){
  p <- ggplot(d, aes_string(x = paste0("`", pah, "`"))) +
    geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
    stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), 
             vjust = -0.5, size = 3.5) +
    labs(title = paste("Distribution of", pah),
         x = paste(pah, "(µg/g)"),
         y = "Frequency") +
    theme_minimal()
  
  print(p)
}

# summarize each PAH
pah_summary_wide <- d %>%
  select(all_of(pah_cols)) %>%
  summarise(across(everything(), 
                   list(total = ~sum(., na.rm = TRUE),
                        n_zero = ~sum(. == 0, na.rm = TRUE),
                        n_nonzero = ~sum(. > 0, na.rm = TRUE)),
                   .names = "{.col}||{.fn}")) %>%
  pivot_longer(everything(),
               names_to = c("PAH", "statistic"),
               names_sep = "\\|\\|") %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  arrange(desc(total))

print(pah_summary_wide)


### Check the LoD reference indicators
# if a measurement is below LoD, indicator='<'
# else indicator='()'
# for example, column 8 (Acenapthene...8) is '<' because Acenapthene (col W) is 0
#   and the LoD is 1.4 (col AM)
# indicators 8:22
# measurements 23:37
# LoDs 39:53
names(d)

# Verify indicator columns (8-22) match PAH measurements (23-37) vs LoDs (39-53)
indicator_cols <- 8:22
pah_measurement_cols <- 23:37
lod_cols <- 39:53

# Create verification table
indicator_check <- tibble(
  column_num = indicator_cols,
  indicator_name = names(d)[indicator_cols],
  pah_name = names(d)[pah_measurement_cols],
  lod_name = names(d)[lod_cols]
)

# Check each PAH
mismatches <- list()

for (i in seq_along(indicator_cols)) {
  indicator_col <- indicator_cols[i]
  pah_col <- pah_measurement_cols[i]
  lod_col <- lod_cols[i]
  
  check_df <- d %>%
    mutate(
      actual_indicator = .[[indicator_col]],
      pah_value = .[[pah_col]],
      lod_value = .[[lod_col]],
      expected_indicator = ifelse(pah_value < lod_value, "<", "()"),
      match = actual_indicator == expected_indicator
    ) %>%
    filter(!match | is.na(match)) %>%
    select(rownum, ParticipantID, SampleID,
           actual_indicator, expected_indicator, pah_value, lod_value)
  
  if (nrow(check_df) > 0) {
    mismatches[[names(d)[pah_col]]] <- check_df
    print(paste("Mismatches found for", names(d)[pah_col]))
    print(check_df)
  }
}

# Summary
if (length(mismatches) == 0) {
  print("All indicator columns match PAH vs LoD comparisons!")
} else {
  print(paste("Mismatches found for", length(mismatches), "PAH(s)"))
}



# 6. Histograms of each PAH LoD
#    Not sure what "below LoD" really means - by each measurement, or some
#    central tendency?

lod_cols <- names(d)[39:53]

for (lod in lod_cols){
  q <- ggplot(d, aes_string(x = paste0("`", lod, "`"))) +
    geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
    stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), 
             vjust = -0.5, size = 3.5) +
    labs(title = paste("Distribution of", lod),
         x = paste(lod, "(µg/g)"),
         y = "Frequency") +
    theme_minimal()
  
  print(q)
}

summary(d[39])
