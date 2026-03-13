# Decision Points & Assumptions for Sensitivity Analysis

This document tracks key analytical decisions and assumptions made during the DFSE Fabric PAH exposure analysis. Each item is a candidate for sensitivity analysis to assess the robustness of conclusions.

---

## Analysis Tracking: Response Variables Across Scripts

| Analysis | Response | Censoring Rate | Script |
|---|---|---|---|
| Simulation validation | Simulated lognormal | ~30% | `00i` |
| Single-PAH MLE vs. naive | Dibenzo(a,h)anthracene | ~89% (400/448 ND) | `01` |
| Baseline model (m1) | `new_totalPAH` (all 15, no imputation) | 86% zeros (384/448) | `02` |
| Improved model (m2) | `totalPAH_imputed` (top 6, MLE imputed) | 0% zeros | `03` |
| QC investigation | Residuals from m2 | — | `03a` |
| Heterogeneous variance model (m3) | `totalPAH_imputed` (top 6, MLE imputed) | 0% zeros | `04` |
| Gamma GLMM, uncollapsed (m4) | `totalPAH_imputed` (top 6, MLE imputed) | 0% zeros — convergence warning | `04` |
| Gamma GLMM, collapsed locations (m4b) | `totalPAH_imputed` (top 6, MLE imputed) | 0% zeros | `04` |

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

**Potential Resolution (adopted for some models on and after m4b):** Collapse ad hoc locations into anatomically adjacent standard locations:

| Ad hoc location | Collapsed into | Rationale |
|---|---|---|
| Shirt (n=1) | Sleeve | PI-approved, same garment region |
| Lower Chest (n=2) | Chest | Adjacent torso locations |
| Fly (n=8) | Pant | Upper leg/crotch region |

This reduces SampleLocation from 13 to 10 levels and eliminates sparse factor levels that caused convergence failures in the Gamma GLMM (m4) and unreliable Cook's distance estimates (decision point #7).

## 7. Cook's distance unreliable with current `SampleLocation` coding

Group-level (participant-level) influence diagnostics via `influence.merMod` are unreliable because removing a participant can eliminate a factor level entirely (e.g., Shirt has n=1). This produces a `rbind` warning about mismatched column lengths and yields a spurious Cook's distance of 25.8 for AA01 — a computational artifact, not a genuine measure of influence. **Resolution required:** collapse `SampleLocation` into broader groups or remove ad hoc locations before running group-level influence diagnostics.

## 8. Heavy right tail in m2 residuals

20 observations have scaled residuals > 2. Investigation (script `03a_qc.R`) shows these are:

- Concentrated at glove locations (Thumb, Finger) and Doffing timing
- Characterized by multiple simultaneous PAH detections (top outliers have 6/6 PAHs detected)
- Distributed across 12 of 23 participants — not concentrated in a few individuals
- Not associated with any particular Burn event

These represent legitimate high-exposure events, not data artifacts. The lognormal residual assumption cannot fully accommodate the tail weight. **Potential actions:** robust mixed models, gamma GLMM with log link, or accept the limitation and note it.

## 9. Reference category for SampleLocation

The current reference level is **Back**, which had zero PAH detections across all analytes — all Back values in the imputed dataset are entirely determined by beta × LOD. This means the model intercept and all location contrasts are relative to imputed noise rather than observed signal. Alternative references to consider: (a) **Chest** — moderate exposure with real detections and reasonable sample size (n=52), (b) a **collapsed body region** group if locations are combined. The choice of reference does not affect model fit or omnibus tests, but changes the interpretation of individual location coefficients.

## 10. Gamma GLMM with location-specific dispersion

Model m4 uses `glmmTMB` with `family = Gamma(link = "log")` and `dispformula = ~SampleLocation` to simultaneously address right-skewed residuals and heterogeneous variance across body locations. The initial fit with all 13 SampleLocation levels produced a convergence warning (`false convergence (8)`) driven by the sparse ad hoc locations (Shirt dispersion SE = 37.3). Collapsing ad hoc locations (decision point #6) is expected to resolve this. **Sensitivity analysis:** compare Gamma GLMM results to the lognormal LME (m3) to assess whether the distributional assumption materially changes the Condition contrasts.

## 11. Location-specific inference is limited by high imputation rates

With ~90% of individual PAH values imputed across most body locations, location-level comparisons are heavily influenced by the imputation assumptions rather than observed data. Locations with real detections (e.g., Neck had 4 outlier observations with genuine PAH signal) may lose statistical significance when their few real values are diluted by the imputed majority. For example, Neck was significantly different from Back in m2 (lognormal LMM, p = 0.002) but non-significant in m4b (Gamma GLMM, p = 0.28) — despite Back having **zero** real detections. The imputed values create an artificial floor that homogenizes locations, making it difficult to distinguish between "truly zero exposure" and "exposed but below LOD." **Implication:** Condition and Timing effects (which cut across all locations) are more robust to this issue than location-specific comparisons.

## 12. Gamma GLMM dispersion estimates invert the variance ranking from the lognormal model

In the lognormal LME (m3), glove locations (Thumb, Finger) had the **highest** residual variance (~1.6–1.8× the reference). In the Gamma GLMM (m4b), these same locations have the **lowest** dispersion (most negative offsets: Thumb = −4.42, Finger = −4.04). This is not contradictory — it reflects different parameterizations. The Gamma separates the mean from the coefficient of variation; glove locations have high means but *consistent* high values, yielding low CV. Meanwhile, Back has the highest dispersion (intercept = 4.12) because its values are entirely imputed noise with no real signal to anchor the mean. **Implication:** Care is needed when comparing variance/dispersion estimates across model families. The substantive question — "which locations have the most variable exposure?" — gets different answers depending on whether variability is measured in absolute terms (lognormal) or relative to the mean (Gamma).

## 13. Visual comparison of raw vs. imputed data via paper doll plots

Create side-by-side paper doll visualizations showing the original `new_totalPAH` (with zeros) alongside `totalPAH_imputed` (MLE beta-substituted). This would communicate the magnitude of the imputation's effect on the data to non-statistical colleagues — how much the zero-dominated landscape changes when non-detects are replaced with small positive values. Adapt plotting code from scripts `00e`/`00f`.

## 14. Hurdle model revisited

Earlier (see discussion in chat history), a hurdle/two-part model was rejected on mechanistic grounds — the zeros reflect measurement sensitivity limits, not a distinct biological process. However, the hurdle model has a pragmatic advantage: **neither part requires imputation**. Part A (logistic: P(detected | covariates)) uses the binary detect/non-detect outcome, which is real data for all 448 observations. Part B (continuous: E[totalPAH | detected, covariates]) models only the ~64 observations with real measured values. This avoids the imputation concerns entirely and preserves the real signals (e.g., glove locations have high exposure, neck sometimes has high exposure, ad hoc fly measurements only in SL/SS conditions). The tradeoff is interpreting two sets of results and the MCAR assumption for Part B (detection is not random — it is associated with location, timing, and condition). Despite this, the hurdle model may tell a richer and more defensible story than the imputation-based approach, particularly when ~90% of values are imputed.


## 15. Duplicate measurements within Participant × Condition × Location × Timing (`Timing_new`)

Script `08` addressed newly identified repeat measurements within the full key (`ParticipantID`, `Condition`, `SampleLocation`, `Timing_new`). Duplicate burden was limited (28 duplicate full-key combinations across 12 participants; 30 extra rows total), but sufficient to justify explicit sensitivity analysis.

Three parallel strategies were fit using the same fixed-effects structure (`Condition + Timing_new + SampleLocation`) with participant random intercepts and location-specific residual variance (`varIdent`):

- **GM collapse:** impute first, then collapse duplicates by geometric mean of `totalPAH_imputed`
- **Median collapse:** impute first, then collapse duplicates by median of `totalPAH_imputed`
- **Repeated-measures random effect:** retain all rows and add nested random effect (`ParticipantID / FullKeyID`)

Findings were highly consistent across all three strategies for key inferential terms (Condition and Timing): effect directions were unchanged, magnitudes were very similar, and significance patterns were stable. The median-collapsed model fit slightly worse than GM collapse; the added duplicate-level random effect did not materially alter fixed-effect inference.

**Proposed reporting frame:** use GM collapse as the parsimonious primary analysis, with median collapse and duplicate-level random effect as prespecified sensitivity analyses demonstrating robustness to duplicate handling.