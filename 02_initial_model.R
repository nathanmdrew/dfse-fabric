### Author: Nathan M. Drew (vom8@cdc.gov)
### Date: 2026-03-05
### Purpose: Begin building a predictive model
### Steps
# 1. Explore the distribution of non-zero new_totalPAH
# 2. Choose c, fit the baseline LMM
# 3. Check residual diagnostics
# 4. Extract variance components and fixed-effect tests

# Required packages (install once):
#install.packages("tidyverse")
#install.packages("lme4")
#install.packages("patchwork")
#install.packages("lmerTest")
#install.packages("emmeans")

library(tidyverse)
library(lme4)
library(patchwork)
library(lmerTest)
library(emmeans)

setwd("C:/Users/vom8/dfse-fabric/")

d <- readRDS(file="data/cleaned_data.RDS")

# histogram and KDE of new_totalPAH
ggplot(d, aes(x=new_totalPAH)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_density(color="red", size=1) +
  labs(title="Histogram and KDE of new_totalPAH", x="new_totalPAH", y="Density") +
  theme_minimal()

# tally zeroes in new_totalPAH
zero_count <- sum(d$new_totalPAH == 0)

table(d$Condition)
table(d$Timing_new)
table(d$`Sample Location`)

# --- Step 2: Choose constant c and fit baseline LMM ---

# Half the minimum non-zero value
c_val <- min(d$new_totalPAH[d$new_totalPAH > 0], na.rm = TRUE) / 2
c_val

# Create log-transformed response
d <- d |>
  mutate(log_totalPAH = log(new_totalPAH + c_val))

# histogram and KDE of new_totalPAH
ggplot(d, aes(x=log_totalPAH)) +
  geom_histogram(aes(y=..density..), bins=30, fill="lightblue", color="black") +
  geom_density(color="red", size=1) +
  labs(title="Histogram and KDE of log_totalPAH", x="log_totalPAH", y="Density") +
  theme_minimal()

# Fit baseline mixed model
m1 <- lmer(
  log_totalPAH ~ Condition + Timing_new + `Sample Location` + (1 | ParticipantID),
  data = d
)

summary(m1)


# --- Step 3: Residual diagnostics for m1 ---

# Extract diagnostics
diag_df <- tibble(
  fitted  = fitted(m1),
  resid   = residuals(m1),
  sresid  = residuals(m1, scaled = TRUE)
)

# 1. Residuals vs. fitted
p1 <- ggplot(diag_df, aes(x = fitted, y = sresid)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(se = FALSE, color = "steelblue", linewidth = 0.8) +
  labs(title = "Residuals vs. Fitted", x = "Fitted values", y = "Scaled residuals") +
  theme_minimal(base_size = 12)

# 2. QQ plot of residuals
p2 <- ggplot(diag_df, aes(sample = sresid)) +
  stat_qq(alpha = 0.4) +
  stat_qq_line(color = "red") +
  labs(title = "Normal Q-Q", x = "Theoretical quantiles", y = "Scaled residuals") +
  theme_minimal(base_size = 12)

# 3. Histogram of residuals
p3 <- ggplot(diag_df, aes(x = sresid)) +
  geom_histogram(bins = 40, fill = "lightblue", color = "black") +
  labs(title = "Residual Distribution", x = "Scaled residuals", y = "Count") +
  theme_minimal(base_size = 12)

# 4. Random effects QQ (participant-level)
re_df <- tibble(
  re = ranef(m1)$ParticipantID[["(Intercept)"]]
)

p4 <- ggplot(re_df, aes(sample = re)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(color = "red") +
  labs(title = "Random Effects Q-Q", x = "Theoretical quantiles",
       y = "Participant intercepts") +
  theme_minimal(base_size = 12)

# Combine
(p1 + p2) / (p3 + p4) +
  plot_annotation(title = "Model m1: Residual Diagnostics")




# --- Step 4: Variance components and fixed-effect tests ---

# Variance components / ICC
vc <- as.data.frame(VarCorr(m1))
icc <- vc$vcov[1] / sum(vc$vcov)

tibble(
  component  = c("ParticipantID", "Residual"),
  variance   = vc$vcov,
  std_dev    = vc$sdcor,
  ICC        = c(icc, NA_real_)
)

# Type III-style F-tests via Satterthwaite (requires lmerTest)
# Refit with lmerTest's lmer to get p-values
m1_lt <- lmerTest::lmer(
  log_totalPAH ~ Condition + Timing_new + `Sample Location` + (1 | ParticipantID),
  data = d
)

anova(m1_lt, type = 3)

# The ANOVA table is striking given how badly the assumptions are violated:

# Sample Location (F = 10.76, p < 2.2e-16) — the dominant signal. This is driven by the glove locations (Finger, Thumb, Palm) having substantially higher values than the body locations. With 12 df in the numerator, there's a lot of heterogeneity across locations.
# Timing_new (F = 71.44, p ≈ 5e-16) — Doffing vs. Donning/Firefighting is the single strongest effect. This makes physical sense — post-fire gear is contaminated.
# Condition (F = 5.80, p = 0.003) — statistically significant, and this is the one most relevant to your hypothesis about SS vs. OL/SL. The SS estimate of 0.43 (t = 3.4) is doing the heavy lifting here; SL is marginal.
# Why this is both encouraging and cautionary:

# The signals are large enough that they're surviving despite the model misspecification. That suggests these are real effects, not artifacts. However, the p-values and F-statistics are unreliable in absolute terms — the bimodal residual structure inflates the residual variance, which could either mask effects (conservative) or distort them (anti-conservative), depending on how the zeros are distributed across factor levels. The direction and relative magnitude of effects are more trustworthy than the exact p-values.

# The ICC of ~8.5% is also informative — it tells you that participant-to-participant differences are modest relative to within-participant variability driven by location, timing, and condition. This is useful context: the exposure story is more about where on the body and when in the firefighting process than about who the firefighter is.

# For the narrative: This baseline model establishes that all three fixed effects are detectable even under a naive modeling approach. The next model iteration (hurdle, restricted to positives, or with proper censoring) should sharpen those estimates and give you defensible inference. The qualitative conclusions are unlikely to change — the effects are too large.

# No, that's not an issue — it's intentional. The ICC (intraclass correlation coefficient) is only meaningful for the grouping factor(s), not the residual. The ICC answers "what proportion of total variance is attributable to between-participant differences?" which is:
# ICC = variance_participant / (variance_participant + variance_residual)
# ICC = 0.085 / (0.085 + 0.92) ≈ 0.085

# The residual variance is the denominator component — it doesn't have its own ICC. The NA just reflects that. If you had multiple random effects (e.g., participant + sample location as random), each would get its own ICC, but the residual row would still be NA.

# Exactly. In plain language: about 8.5% of the total variation in log-transformed total PAH exposure is due to differences between firefighters, while the remaining ~91.5% is due to differences within the same firefighter — i.e., variation driven by where on the body the sample was taken, when it was collected, and the experimental condition.

# Put another way: knowing which firefighter a measurement came from tells you relatively little about the expected PAH level. The context of the measurement (location, timing, condition) matters far more than who the person is.

#----------- Step 5: Post-hoc contrasts for Condition -----------
d <- d |>
  mutate(`Sample Location` = factor(`Sample Location`))

# Estimated marginal means for Condition
emm_cond <- emmeans(m1_lt, "Condition")
emm_cond

# Pairwise contrasts (SS vs OL, SS vs SL, SL vs OL)
contrast(emm_cond, method = "pairwise", adjust = "tukey")

#---------- Step 6: Save model and results for next steps -----------
# Model objects
saveRDS(m1, "02_output/m1_baseline_lmer.rds")
saveRDS(m1_lt, "02_output/m1_baseline_lmerTest.rds")

# Diagnostics
saveRDS(diag_df, "02_output/m1_diag_df.rds")
saveRDS(re_df, "02_output/m1_re_df.rds")

# Variance components and ICC
saveRDS(vc, "02_output/m1_variance_components.rds")
saveRDS(icc, "02_output/m1_icc.rds")

# Contrasts
saveRDS(emm_cond, "02_output/m1_emmeans_condition.rds")

# Scalar used for log transform
saveRDS(c_val, "02_output/c_val.rds")

# Diagnostic plots
ggsave("02_output/m1_diagnostic_plots.png",
       (p1 + p2) / (p3 + p4) + plot_annotation(title = "Model m1: Residual Diagnostics"),
       width = 12, height = 8, dpi = 300)
