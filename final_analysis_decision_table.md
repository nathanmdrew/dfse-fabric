# Final Analysis Decision Table

## Primary analysis specification

| Domain | Decision | Rationale | Notes for manuscript |
|---|---|---|---|
| Outcome construction | `totalPAH_imputed` from **Top 4 PAHs** (global detect-rate ranking) | Balances signal coverage with lower imputation burden than Top 6/15 | Top 1, Top 6, Top 15 retained as sensitivity/boundary runs |
| Non-detect handling | **PAH-level MLE beta** via `EnvStats::elnormCensored`, impute ND as `beta_j * LOD_i` | Consistent with prior scripts/manuscript logic; interpretable | Beta estimated once per PAH; reference `LOD` = median non-missing LOD for that PAH |
| Repeated measures | Collapse replicate rows by geometric mean within `ParticipantID × Condition × Timing_new × SampleLocation` | Preserves within-cell central tendency; aligns with prior scripts | Collapse performed after imputation at PAH level |
| Model form | `nlme::lme(log_totalPAH_imp ~ Condition + Timing_new + SampleLocation, random = ~1 | ParticipantID, weights = varIdent(~1|SampleLocation), method = "REML")` | Matches established model structure and accommodates heteroscedasticity by location | Additive form only (no interactions) due to sample size constraints |
| Primary estimand | Pairwise **Condition** contrasts (`OL-SL`, `OL-SS`, `SL-SS`) from `emmeans`, Tukey-adjusted | Directly addresses PI hypothesis about SS | Timing contrast treated as secondary |
| Secondary estimand | Pairwise `Timing_new` contrast from `emmeans` | Complements Condition results | Adjustment rule pre-specified and reported |

## Sensitivity / boundary analyses

| Lever | Levels evaluated | Purpose | Reporting format |
|---|---|---|---|
| PAH set size | `top1` through `top10`, plus `top15` | Boundary assessment for response definition | Direction stability, estimate range, significance frequency |
| Imputation factor | Beta sensitivity previously evaluated (script 09) | Quantify dependence on ND substitution magnitude | Reported separately from PAH-count sensitivity |
| Non-imputed baseline | As-is and detects-only variants (script 10) | Lower-bound / alternate handling of NDs | Tiered comparability labels retained |

## Comparability tiers for interpretation

| Tier | Definition | Use in conclusions |
|---|---|---|
| Tier 1 (full) | Heteroscedastic `varIdent` model fit succeeds | Primary cross-run comparisons |
| Tier 2 (fallback) | Homoscedastic fallback required | Sensitivity context only |
| Tier 3 | Model fails | Not used for inferential comparison |

## Pre-specified reviewer-facing statements

1. Direction/magnitude stability is emphasized over single-run p-values due to limited signal.
2. Imputation uncertainty is **not** propagated into standard errors; inference is conditional on imputed values.
3. Primary narrative is based on Condition contrasts under the pre-specified primary model; all alternatives are sensitivity analyses.
4. PAH ranking is fixed globally from raw detect rates before modeling.

## Final analysis run plan (clean script)

1. Load data and harmonize factors/location recode.
2. Build fixed PAH ranking and define sets.
3. Run primary model (Top 4, MLE-beta imputation, GM collapse, heteroscedastic LME).
4. Extract primary and secondary contrasts.
5. Run sensitivity matrix (Top 1/4/6/15 minimum; optionally Top 1–10 + 15).
6. Produce boundary summary tables and plots.
7. Save `.rds` artifacts and `.png` figures with manifest.
