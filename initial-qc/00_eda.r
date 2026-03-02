library(dplyr)
library(readxl)
library(tidyr)
library(stringr)

d <- read_excel("Human_Subject_Results_Ver5_master.xlsx", sheet="Fabric")
d <- as.data.frame(d)

# Parse SampleID into components (keep SampleID as the original full value)
d <- d %>%
  mutate(SampleID = as.character(SampleID)) %>%
  mutate(ParticipantID = str_extract(SampleID, "^[^_]+")) %>%
  mutate(rest = str_remove(SampleID, "^[^_]+_")) %>%
  separate(rest, into = c("Condition2", "SampleLocation2", "SampleID2_raw"),
           sep = "_", remove = TRUE, extra = "merge", fill = "right") %>%
  mutate(SampleID2 = as.integer(str_extract(SampleID2_raw, "\\d+"))) %>%
  select(-SampleID2_raw)

# Decode SampleLocation2
d <- d %>%
  mutate(
    SampleLocation_decoded = case_when(
      SampleLocation2 == "C" ~ "Chest",
      SampleLocation2 == "GF" ~ "Finger",
      SampleLocation2 == "GP" ~ "Palm",
      SampleLocation2 == "GT" ~ "Thumb",
      SampleLocation2 == "LP" ~ "Lower Pant",
      SampleLocation2 == "N" ~ "Neck",
      SampleLocation2 == "P" ~ "Pant",
      SampleLocation2 == "SO" ~ "Sock",
      SampleLocation2 == "S" ~ "Shirt",  # Note: "S" may correspond to Shirt or Sleeve
      SampleLocation2 == "B" ~ "Back",
      SampleLocation2 == "F" ~ "Fly",
      SampleLocation2 == "LC" ~ "Lower Chest",
      SampleLocation2 == "GL" ~ "Finger",
      TRUE ~ NA_character_
    )
  )

# Create match flags
d <- d %>%
  mutate(
    Condition_match = ifelse(!is.na(Condition) & !is.na(Condition2) &
                             as.character(Condition) == as.character(Condition2), TRUE, FALSE),
    SampleLocation_match = ifelse(
      !is.na(`Sample Location`) & !is.na(SampleLocation2),
      ifelse(SampleLocation2 == "S",
             `Sample Location` %in% c("Shirt", "Sleeve"),
             as.character(`Sample Location`) == SampleLocation_decoded),
      FALSE
    )
  )

# Summary of mismatches (Condition or SampleLocation)
mismatches <- d %>% filter(!(Condition_match & SampleLocation_match))
print(paste("Total rows:", nrow(d)))
print(paste("Mismatching rows (Condition or SampleLocation):", nrow(mismatches)))
if(nrow(mismatches) > 0) {
  print(head(mismatches %>% select(SampleID, Condition, Condition2, `Sample Location`, SampleLocation2, SampleLocation_decoded, ParticipantID, SampleID2), 50))
}

#print(names(d))

#print(unique(d$SampleLocation2))

# Save d as RDS
saveRDS(d, "fabric_data.rds")

# Export d as CSV
#write.csv(d, "fabric_data.csv", row.names = FALSE)

print(unique(d$ParticipantID))

qc <- unique(d$ParticipantID)
qc2 <- str_sub(qc, -2)

print(unique(qc2))
