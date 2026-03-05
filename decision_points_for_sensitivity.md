# Decision Points & Assumptions for Sensitivity Analysis

This document tracks key analytical decisions and assumptions made during the DFSE Fabric PAH exposure analysis. Each item is a candidate for sensitivity analysis to assess the robustness of conclusions.

---

## 1. Top 6 PAHs used to compute the sum

The composite exposure variable (`totalPAH_imputed`) is computed as the sum of the top 6 PAHs ranked by detection count:

- Dibenzo(a,h)anthracene
- Benzo(b)fluoranthene
- Indeno(1,2,3-cd)pyrene
- Benzo(a)pyrene
- Benzo(g,h,i)perylene
- Benzo(a)anthracene

These were verified against the previous analysis (Wilkinson et al. 2025). The remaining 9 PAHs are excluded. **Sensitivity analysis:** include all 15 PAHs with imputation applied to each.

## 2. Median LOD as the representative LOD

A single representative LOD per PAH (the median of all sample-level LODs) is used in the MLE beta derivation. Individual samples have varying LODs. **Sensitivity analysis:** use per-sample LOD values instead of a single median.

## 3. MLE censored lognormal assumed for each PAH

The `elnormCensored()` function from `EnvStats` is used with `method = "mle"`, assuming a lognormal distribution for each PAH. **Sensitivity analysis:** explore alternative distributional assumptions (gamma, Weibull, etc.).

## 4. ~~Constant `c` in log(y + c)~~ — Not needed after imputation

The baseline model (m1) required a constant `c = min(nonzero)/2` for the log transform due to 384/448 zeros in `new_totalPAH`. After MLE beta-substitution imputation, all values of `totalPAH_imputed` are positive, so `log()` can be applied directly. This decision point is resolved.

## 5. Top 6 PAHs only in `totalPAH_imputed`

`totalPAH_imputed` sums only the top 6 PAHs, while the original `new_totalPAH` summed all 15. The top 6 capture ~90% of the total signal at the maximum, but some location-specific signals (e.g., Lower Pant) may come from PAHs outside the top 6. **Sensitivity analysis:** recompute with all 15 PAHs imputed.

## 6. Ad hoc sample locations (Shirt, Lower Chest, Fly)

Shirt (n=1), Lower Chest (n=2), and Fly (n=8) were collected ad hoc due to visible contamination at sample collection time — different criteria than the standard protocol locations. Shirt could likely be combined with Sleeve (per the project PI). **Sensitivity analysis:** (a) combine Shirt with Sleeve, (b) omit all ad hoc locations, (c) otherwise account for the design difference.

## 7. Cook's distance unreliable with current `SampleLocation` coding

Group-level (participant-level) influence diagnostics via `influence.merMod` are unreliable because removing a participant can eliminate a factor level entirely (e.g., Shirt has n=1). This produces a `rbind` warning about mismatched column lengths and yields a spurious Cook's distance of 25.8 for AA01 — a computational artifact, not a genuine measure of influence. **Resolution required:** collapse `SampleLocation` into broader groups or remove ad hoc locations before running group-level influence diagnostics.

## 8. Heavy right tail in m2 residuals

20 observations have scaled residuals > 2. Investigation (script `03a_qc.R`) shows these are:

- Concentrated at glove locations (Thumb, Finger) and Doffing timing
- Characterized by multiple simultaneous PAH detections (top outliers have 6/6 PAHs detected)
- Distributed across 12 of 23 participants — not concentrated in a few individuals
- Not associated with any particular Burn event

These represent legitimate high-exposure events, not data artifacts. The lognormal residual assumption cannot fully accommodate the tail weight. **Potential actions:** robust mixed models, gamma GLMM with log link, or accept the limitation and note it.