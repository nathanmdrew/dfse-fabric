# Plan for Handling Left-Censored PAH Data (Below Limit of Detection)

## Overview
The Total_PAH measurements contain values recorded as 0 when the sum of individual PAH compounds falls below the detection limit. Individual PAH columns use "<" to indicate values below their specific LOD, with LOD values stored in corresponding `_PAH_LOD` columns. This creates a hierarchical censoring structure that requires careful handling.

## Phase 1: Data Documentation & Structure Identification

### Objective
Map individual PAH measurement columns to their corresponding LOD columns and flag censored observations.

### Steps
```r
# 1. Identify all PAH measurement and LOD columns
pah_cols <- names(d) %>% 
  str_subset("^[A-Z].*[0-9]$")  # Adjust regex to match your PAH column naming

lod_cols <- names(d) %>% 
  str_subset("_LOD$")

# 2. Examine the structure
print(pah_cols)
print(lod_cols)

# 3. Verify mapping between PAH columns and LOD columns
# For each PAH column, identify its corresponding LOD column
```

Output
- Complete list of individual PAH measurement columns
- Complete list of LOD columns
- Mapping document showing which PAH goes with which LOD

---

## Phase 2: Flag Individual Censored Values
Objective
Identify which individual PAH measurements are below their respective LODs.

Steps
```r
# 1. Flag rows where PAH value is "<" (censored)
d <- d %>%
  mutate(across(all_of(pah_cols), 
                list(is_censored = ~str_detect(., "<")),
                .names = "{.col}_is_censored"))

# 2. Convert PAH values to numeric, preserving NA for censored values
d <- d %>%
  mutate(across(all_of(pah_cols),
                ~ifelse(str_detect(., "<"), NA_real_, as.numeric(.)),
                .names = "{.col}_numeric"))

# 3. For censored values, retrieve the LOD
# Example structure (adjust based on actual column names):
d <- d %>%
  mutate(
    Acenaphthene_imputed = case_when(
      Acenaphthene_is_censored ~ Acenaphthene_LOD,
      !is.na(Acenaphthene_numeric) ~ Acenaphthene_numeric,
      TRUE ~ NA_real_
    )
  )
```

Output
- Flag columns indicating which individual PAH values are censored
- Numeric versions of PAH columns for analysis
- Decision on whether to use LOD or LOD/2 for censored values

---

## Phase 3: Handle Total_PAH Censoring
Objective
Properly characterize Total_PAH as censored when any individual component is censored.

Key Concept
When summing left-censored values, the total becomes left-censored at the sum of individual LODs. A Total_PAH of 0 actually means: Total_PAH < (sum of all component LODs).

Steps
```r
# 1. Calculate which rows have any censored components
d <- d %>%
  mutate(
    n_components_censored = rowSums(select(., ends_with("_is_censored"))),
    any_component_censored = n_components_censored > 0,
    sum_of_lods = rowSums(select(., all_of(lod_cols)), na.rm = TRUE)
  )

# 2. Recalculate Total_PAH from imputed component values
# (Using LOD or LOD/2 for censored values)
d <- d %>%
  mutate(
    total_pah_calculated = rowSums(select(., ends_with("_imputed")), na.rm = TRUE)
  )

# 3. Flag censoring status of total
d <- d %>%
  mutate(
    total_pah_status = case_when(
      is.na(total_pah_calculated) ~ "Missing",
      any_component_censored ~ paste0("Censored (< ", round(sum_of_lods, 2), ")"),
      total_pah_calculated == 0 ~ "All components < LOD",
      TRUE ~ "Detected"
    ),
    total_pah_is_censored = any_component_censored | total_pah_calculated == 0
  )

# 4. Summary of censoring
table(d$total_pah_status)
```

Output
- Flag indicating which Total_PAH values are censored
- Censoring limits (sum of LODs) for each observation
- Summary table of censoring patterns

---

## Phase 4: Choose Analysis Approach (Future Implementation)
Option 1: Maximum Likelihood (Recommended for inference)
- Packages: censReg, survival, or Bayesian approaches
- Advantages: Treats zeros as left-censored without imputation; preserves statistical properties; appropriate for dose-response analysis
- Best for: Regression modeling, BMD analysis, hypothesis testing

Option 2: Multiple Imputation
- Packages: mice, mi
- Method: Impute censored values from a truncated normal/lognormal distribution (lower bound = 0, upper bound = sum of LODs)
- Advantages: Allows use of standard statistical methods
- Disadvantages: More computationally intensive; loss of information about censoring

Option 3: Substitution Methods (Not recommended for inference)
- Methods: LOD/2 or LOD/√2 replacement
- Advantages: Simple, no special software needed
- Disadvantages: Biased estimates; inappropriate for statistical inference
Use only for: Descriptive statistics

---

## Next Steps
- Immediately: Complete Phase 1—identify all PAH and LOD columns and verify the mapping
- In initial EDA: Implement Phase 2 and 3 to characterize censoring patterns by subject and location
- Before modeling: Evaluate which analysis approach (ML, MI, or other) best suits your research questions and available software
- Documentation: Keep track of imputation decisions (LOD vs. LOD/2) and justify the choice

---

## Questions to Resolve
- Which censoring substitution should be used if needed? (LOD or LOD/2 or other?)
- Are there any missing values in the LOD columns themselves?
- Do censoring patterns differ significantly by sample location or subject?
- What is the expected proportion of censored observations in the final dataset?