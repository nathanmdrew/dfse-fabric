## Summary of Ganser & Hewett (2010)

**"An Accurate Substitution Method for Analyzing Censored Data"**
*Journal of Occupational and Environmental Hygiene, 7(4): 233–244*

---

### Problem Being Addressed

In occupational and environmental hygiene, exposure datasets frequently contain **left-censored observations** — measurements reported only as "less than the limit of detection" (LOD). Accurately estimating distributional parameters (GM, GSD, 95th percentile, and Mean) from such datasets is a core analytical challenge.

---

### Methods Compared

The paper evaluates four approaches for handling censored data:

1. **LOD/2 substitution** — replace each non-detect with LOD × 0.5
2. **LOD/√2 substitution** — replace each non-detect with LOD × 0.707
3. **Maximum Likelihood Estimation (MLE)** — considered the "gold standard," uses the Newton-Raphson method to maximize a likelihood function over both detected and censored observations
4. **β-substitution (new method)** — a moment-based substitution method that calculates a dataset-specific adjustment factor (β) from the uncensored observations to replace each LOD value

---

### The β-Substitution Algorithm (Key Innovation)

Rather than using a fixed factor (1/2 or 1/√2), β-substitution computes an **empirically derived β factor** separately for estimating the **Mean** (β_MEAN) and **GM** (β_GM). The algorithm proceeds as follows:

1. Compute the mean of log-transformed detects (ȳ) and the Z-value corresponding to the censoring fraction (k/n)
2. Estimate ŝy (initial estimate of ln(GSD)) from the uncensored data
3. Calculate β_MEAN and β_GM using closed-form expressions involving normal PDF/CDF functions
4. Substitute non-detects with β·LOD and compute sample statistics
5. Derive GSD from the ratio of mean to GM, then compute the 95th percentile

For **multiple LODs**, the method uses a geometric mean of the individual LODs, weighted by the number of non-detects at each LOD. The algorithm is fully implementable in a spreadsheet (Appendix B provides Excel cell formulas).

---

### Simulation Study Design

Performance was evaluated across **100,000 simulated datasets** per condition, spanning:

| Parameter | Range |
|---|---|
| Sample size (Sim 1) | 20–100 |
| Sample size (Sim 2) | 5–19 |
| True GSD | 1.2–4.0 |
| Percent censored | 1%–50% |
| Distribution type | Single lognormal; contaminated (mixture) lognormal |
| Number of LODs | Single and multiple |

Performance metrics were **bias** and **root mean square error (rMSE)**.

---

### Key Findings

- **β-substitution closely tracks MLE** in both bias and rMSE across all simulation scenarios — for GM, GSD, 95th percentile, and Mean
- **For small samples (n = 5–19)**, β-substitution was *superior* to MLE in terms of both bias and rMSE, particularly for the 95th percentile
- **LOD/2 and LOD/√2** showed highly variable bias — strongly positive or negative depending on true GSD — and this bias did not diminish with increasing sample size
- β-substitution handles **multiple LODs and non-lognormal (contaminated) distributions** robustly
- At high censoring (50–80%), β-substitution continued to perform comparably to MLE

---

### Relevance to Your R Coding Context

For a NIOSH biostatistician working with occupational exposure data (including nanoscale particle data), this paper is directly relevant when:

- Datasets contain **non-detects or values below an analytical LOD**
- You need to estimate **lognormal distribution parameters** (GM, GSD), **upper percentiles** (e.g., 95th), or the **arithmetic mean** from censored data
- MLE is unavailable or impractical (e.g., large databases, job-exposure matrices)
- Sample sizes are **small (n < 20)**, where MLE can be severely biased

The β-substitution algorithm can be implemented in R using standard `stats` functions (`pnorm`, `dnorm`, `qnorm`) and is well-suited for integration into a `tidyverse`-based analysis pipeline.