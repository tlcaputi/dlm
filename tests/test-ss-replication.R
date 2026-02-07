#!/usr/bin/env Rscript
# test-ss-replication.R
#
# Verify the dlm package produces results IDENTICAL to:
# 1. fixest event study with binned endpoints (on the same sample)
# 2. The Schmidheiny & Siegloch (2023, JAE) replication using their bf2017 data
#
# This is the gold standard test: machine epsilon equivalence.

suppressPackageStartupMessages({
  library(haven)
  library(data.table)
  library(dplyr)
  library(fixest)
  library(plm)
  library(lmtest)
  library(zoo)
  library(glue)
  library(rlang)
  library(scales)
  library(ggplot2)
  library(tidyr)
  library(logger)
})

# Source the dlm package functions directly
for (f in list.files("/tmp/dlm-work/R", pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

cat("\n========================================\n")
cat("TEST 1: DLM-ES Equivalence on generate_data()\n")
cat("========================================\n\n")

# Already proven in previous session, but let's reconfirm
df <- generate_data(seed = 42)
from_rt <- -3; to_rt <- 3; ref_period <- -1

dlm_results <- distributed_lags_models(
  data = df, exposure_data = df,
  from_rt = from_rt, to_rt = to_rt,
  outcomes = "outcome", exposure = "post",
  unit = "group", time = "time",
  ref_period = ref_period
)
dlm_betas <- dlm_results[[1]]$betas
dlm_N <- nobs(dlm_results[[1]]$model)

# Restrict ES to same sample (critical!)
tw_min <- to_rt + 1
tw_max <- max(df$time) - (abs(from_rt) - 1)
df_r <- df[df$time >= tw_min & df$time <= tw_max, ]
stopifnot(nrow(df_r) == dlm_N)

es_model <- feols(
  outcome ~ i(years_to_treatment, treat,
              ref = c(ref_period, -1000),
              bin = list("-3+" = ~x <= -3, "3+" = ~x >= 3)) | group + time,
  data = df_r, cluster = ~group
)
es_ct <- as.data.frame(es_model$coeftable)

max_coef_diff <- max(abs(dlm_betas$coef - es_ct[, 1]))
max_se_diff <- max(abs(dlm_betas$se - es_ct[, 2]))

cat("  Max coefficient difference:", format(max_coef_diff, scientific = TRUE), "\n")
cat("  Max SE difference:         ", format(max_se_diff, scientific = TRUE), "\n")
cat("  PASS:", max_coef_diff < 1e-10 && max_se_diff < 1e-10, "\n")

cat("\n========================================\n")
cat("TEST 2: DLM-ES Equivalence on BF2017 Data (S&S Case 2)\n")
cat("========================================\n\n")

# Load the Baker & Fradkin (2017) data used by S&S
bf2017 <- read_dta("/tmp/ss-replication/resources/bf2017.dta")
bf2017 <- as.data.frame(bf2017)

# Create a SEQUENTIAL time index (critical for our package)
# yearmonth like 200601 has gaps (200612 -> 200701) which breaks tidyr::complete()
bf2017 <- bf2017[order(bf2017$state, bf2017$year, bf2017$month), ]
time_periods <- sort(unique(paste0(bf2017$year, sprintf("%02d", bf2017$month))))
time_map <- data.frame(
  ym_str = time_periods,
  time_idx = seq_along(time_periods),
  stringsAsFactors = FALSE
)
bf2017$ym_str <- paste0(bf2017$year, sprintf("%02d", bf2017$month))
bf2017 <- merge(bf2017, time_map, by = "ym_str")
bf2017 <- bf2017[order(bf2017$state, bf2017$time_idx), ]

# Create treatment adoption indicator (increase by >= 13 weeks)
bf2017 <- bf2017 %>%
  group_by(state) %>%
  mutate(
    D_PBD_incr = ifelse(PBD - dplyr::lag(PBD) >= 13, 1, 0),
    D_PBD_incr = ifelse(is.na(D_PBD_incr), 0, D_PBD_incr)
  ) %>%
  ungroup()

# Create treatment status variable (cumulative sum of adoption indicator)
bf2017 <- bf2017 %>%
  group_by(state) %>%
  mutate(PBD_incr = cumsum(D_PBD_incr)) %>%
  ungroup()

# Create covariate
bf2017$frac_total_ui <- (bf2017$cont_claims + bf2017$init_claims) / bf2017$population

# Crisis sample: year <= 2011
bf2017_crisis <- bf2017[bf2017$year <= 2011, ]

# Create log outcome
bf2017_crisis$log_GJSI <- log(bf2017_crisis$GJSI)

cat("BF2017 crisis sample: ", nrow(bf2017_crisis), " obs\n")
cat("States: ", length(unique(bf2017_crisis$state)), "\n")
cat("Yearmonths: ", length(unique(bf2017_crisis$yearmonth)), "\n\n")

# --- Run DLM using our package ---
cat("Running DLM via our package...\n")

# Our package needs: unit, time, exposure (the LEVEL variable)
# from_rt = -3 means 2 leads (lead2, lead1) + contemporaneous + 4 lags = from_rt=-3, to_rt=4
# Wait - S&S uses effect window -3 to +4, which means:
#   leads: period -3, -2 (2 leads of the level variable)
#   ref: period -1 (zero by construction)
#   lags: period 0, 1, 2, 3, 4 (contemporaneous + 4 lags)
# In our package: from_rt = -3, to_rt = 4
# This creates leads: abs(from_rt)-1 = 2 leads (lead2, lead1)
# And lags: 0 to to_rt = lag0, lag1, lag2, lag3, lag4
# Total gamma params: 2 + 5 = 7 (matching S&S)
# Before periods (< ref_period=-1): from_rt to ref_period-1 = -3, -2 → 2 periods
# After periods (> ref_period=-1): 0, 1, 2, 3, 4 → 5 periods

dlm_bf <- distributed_lags_model(
  data = bf2017_crisis,
  exposure_data = bf2017_crisis,
  from_rt = -3, to_rt = 4,
  outcome = "log_GJSI",
  exposure = "PBD_incr",
  unit = "state", time = "time_idx",
  covariates = c("frac_total_ui"),
  ref_period = -1
)

dlm_bf_betas <- dlm_bf$betas
dlm_bf_N <- nobs(dlm_bf$model)
dlm_bf_gamma <- dlm_bf$model$coefficients[1:7]
dlm_bf_vcov <- dlm_bf$vcov

cat("\nDLM gamma coefficients:\n")
print(dlm_bf_gamma)
cat("\nDLM beta coefficients (cumsum-transformed):\n")
print(dlm_bf_betas)
cat("DLM N:", dlm_bf_N, "\n\n")

# --- Run S&S-style DLM using plm (for comparison) ---
cat("Running S&S-style DLM via plm...\n")

# Prepare plm panel using yearmon (as S&S do)
bf2017_plm <- bf2017
bf2017_plm$yearmonth_ym <- as.yearmon(paste(bf2017_plm$year, sprintf("%02d", bf2017_plm$month)), "%Y%m")
bf2017_plm <- pdata.frame(bf2017_plm, index = c("state", "yearmonth_ym"))

# Estimate DLM in levels (exact S&S specification)
estim_ss <- plm(log(GJSI) ~
                  lead(PBD_incr, k=2) + lead(PBD_incr, k=1) +
                  PBD_incr +
                  lag(PBD_incr, k=1) + lag(PBD_incr, k=2) + lag(PBD_incr, k=3) + lag(PBD_incr, k=4) +
                  frac_total_ui + as.factor(yearmonth_ym),
                data = bf2017_plm,
                subset = (year <= 2011),
                effect = "individual",
                model = "within")

gamma_ss <- estim_ss$coefficients[1:7]
vcov_ss <- vcovHC(estim_ss, type = "sss", cluster = "group")[1:7, 1:7]

# S&S cumulative sum transformation
beta_ss_coef <- c(-revcumsum(gamma_ss[1:2]), 0, cumsum(gamma_ss[3:7]))
beta_ss_se <- c(abs(serevcumsum(vcov_ss[1:2, 1:2])), 0, secumsum(vcov_ss[3:7, 3:7]))

cat("\nS&S gamma coefficients (plm):\n")
print(gamma_ss)
cat("\nS&S beta coefficients (plm):\n")
print(data.frame(time = -3:4, coef = beta_ss_coef, se = beta_ss_se))
cat("S&S N:", nobs(estim_ss), "\n\n")

# --- Compare DLM (our package) vs S&S DLM (plm) ---
cat("--- Comparison: Our DLM vs S&S DLM (plm) ---\n")
cat("(Different packages, so small SE differences expected due to df corrections)\n\n")

# Point estimates should be very close
gamma_diff <- max(abs(dlm_bf_gamma - gamma_ss))
cat("Max gamma difference:", format(gamma_diff, scientific = TRUE), "\n")

# Our betas vs S&S betas (excluding ref period = 0)
our_coefs <- c(dlm_bf_betas$coef[1:2], dlm_bf_betas$coef[3:7])
ss_coefs <- c(beta_ss_coef[1:2], beta_ss_coef[4:8])
coef_diff <- max(abs(our_coefs - ss_coefs))
cat("Max beta coefficient difference:", format(coef_diff, scientific = TRUE), "\n")

# If N differs, note it
if (dlm_bf_N != nobs(estim_ss)) {
  cat("WARNING: Sample sizes differ! DLM:", dlm_bf_N, "vs S&S:", nobs(estim_ss), "\n")
  cat("This is expected if fixest and plm handle boundary obs differently.\n")
} else {
  cat("Sample sizes match:", dlm_bf_N, "\n")
}

cat("\n========================================\n")
cat("TEST 3: DLM vs ES on BF2017 (fixest only, same-sample)\n")
cat("  This is the TRUE equivalence test\n")
cat("========================================\n\n")

# For the fixest ES, we need treatment adoption indicators with binned endpoints
# following S&S eq. (5)

# First, let's get the exact sample used by our DLM
# Our DLM drops obs where leads/lags are NA
# For from_rt=-3, to_rt=4: need 2 leads and 4 lags
# So effective sample is: time >= min_time + 4, time <= max_time - 2

# Get the DLM's actual sample size to match
cat("DLM sample size:", dlm_bf_N, "\n")

# Create treatment adoption indicator leads/lags for the ES
# The ES uses the adoption indicators (first differences), not levels
bf_es <- bf2017_crisis %>%
  group_by(state) %>%
  arrange(time_idx) %>%
  mutate(
    # Leads of adoption indicator
    F3_D = dplyr::lead(D_PBD_incr, 3),
    F2_D = dplyr::lead(D_PBD_incr, 2),
    F1_D = dplyr::lead(D_PBD_incr, 1),
    # Lags
    L1_D = dplyr::lag(D_PBD_incr, 1),
    L2_D = dplyr::lag(D_PBD_incr, 2),
    L3_D = dplyr::lag(D_PBD_incr, 3),
    L4_D = dplyr::lag(D_PBD_incr, 4)
  ) %>%
  ungroup()

# Determine the DLM sample by checking which obs have non-NA leads/lags
# The DLM uses 2 leads of PBD_incr and 4 lags
bf_dlm_check <- bf2017_crisis %>%
  group_by(state) %>%
  arrange(time_idx) %>%
  mutate(
    lead2 = dplyr::lead(PBD_incr, 2),
    lag4 = dplyr::lag(PBD_incr, 4)
  ) %>%
  ungroup() %>%
  filter(!is.na(lead2) & !is.na(lag4))

cat("DLM effective sample (calculated):", nrow(bf_dlm_check), "\n")

# Restrict ES to same sample, then build binned endpoints
# Following S&S eq. (5)
bf_es_r <- bf_es %>%
  semi_join(bf_dlm_check, by = c("state", "time_idx")) %>%
  group_by(state) %>%
  arrange(time_idx) %>%
  mutate(
    # Binned far-lag endpoint: cumsum of L4
    L4bin_D = cumsum(ifelse(is.na(L4_D), 0, L4_D)),
    # Binned far-lead endpoint: reverse cumsum of F3
    F3bin_D = revcumsum(ifelse(is.na(F3_D), 0, F3_D))
  ) %>%
  ungroup()

cat("ES restricted sample:", nrow(bf_es_r), "\n")

# Estimate ES with fixest using binned endpoints
# Regressors: F3bin, F2, [ref=-1 omitted], D_PBD_incr, L1, L2, L3, L4bin
# Plus covariate + FEs
es_bf <- feols(
  log_GJSI ~ F3bin_D + F2_D + D_PBD_incr + L1_D + L2_D + L3_D + L4bin_D + frac_total_ui | state + time_idx,
  data = bf_es_r, cluster = ~state
)

es_bf_ct <- as.data.frame(es_bf$coeftable)
es_bf_coef <- es_bf_ct[1:7, 1]  # 7 ES coefficients
es_bf_se <- es_bf_ct[1:7, 2]

cat("\nES coefficients (fixest):\n")
print(data.frame(
  period = c(-3, -2, 0, 1, 2, 3, 4),
  coef = es_bf_coef,
  se = es_bf_se
))
cat("ES N:", nobs(es_bf), "\n\n")

# Now compare: DLM betas should equal ES betas
# DLM betas: periods -3, -2, 0, 1, 2, 3, 4 (ref period -1 excluded)
cat("DLM beta coefficients:\n")
print(dlm_bf_betas)

# Map: DLM period -3 ↔ ES F3bin, -2 ↔ F2, 0 ↔ D_PBD_incr, 1 ↔ L1, etc.
dlm_coefs_ordered <- dlm_bf_betas$coef
dlm_ses_ordered <- dlm_bf_betas$se

cat("\n--- EQUIVALENCE CHECK ---\n")
cat("DLM coefs:  ", paste(format(dlm_coefs_ordered, digits = 10), collapse = ", "), "\n")
cat("ES  coefs:  ", paste(format(es_bf_coef, digits = 10), collapse = ", "), "\n")

coef_diffs <- abs(dlm_coefs_ordered - es_bf_coef)
se_diffs <- abs(dlm_ses_ordered - es_bf_se)

cat("\nCoefficient differences:\n")
print(coef_diffs)
cat("Max coef diff:", format(max(coef_diffs), scientific = TRUE), "\n")

cat("\nSE differences:\n")
print(se_diffs)
cat("Max SE diff:", format(max(se_diffs), scientific = TRUE), "\n")

if (max(coef_diffs) < 1e-10 && max(se_diffs) < 1e-10) {
  cat("\n*** PASS: DLM = ES to machine precision on BF2017 data ***\n")
} else if (max(coef_diffs) < 1e-10) {
  cat("\n*** PARTIAL: Coefficients match but SEs differ ***\n")
  cat("(This may be due to sample size differences between DLM and ES)\n")
} else {
  cat("\n*** INVESTIGATING: Results don't match to machine precision ***\n")
  cat("N DLM:", dlm_bf_N, "  N ES:", nobs(es_bf), "\n")
}

cat("\n========================================\n")
cat("TEST 4: Verify secumsum/serevcumsum match S&S functions exactly\n")
cat("========================================\n\n")

# Create a test vcov matrix
set.seed(42)
A <- matrix(rnorm(25), 5, 5)
test_vcov <- A %*% t(A)  # positive definite

# Our secumsum
our_secumsum <- secumsum(test_vcov)

# S&S secumsum (identical algorithm)
ss_secumsum <- function(vcov) {
  L <- dim(vcov)[1]
  se <- c()
  for (i in c(1:L)) {
    a <- matrix(rep(1, i), nrow = 1)
    V <- a %*% vcov[1:i, 1:i] %*% t(a)
    se[i] <- sqrt(V)
  }
  return(se)
}

ss_result <- ss_secumsum(test_vcov)
diff_secumsum <- max(abs(our_secumsum - ss_result))
cat("secumsum max diff:", format(diff_secumsum, scientific = TRUE), "\n")

# Our serevcumsum
our_serevcumsum <- serevcumsum(test_vcov)

# S&S serevcumsum
ss_serevcumsum <- function(vcov) {
  L <- dim(vcov)[1]
  se <- c()
  for (i in c(L:1)) {
    a <- matrix(rep(1, L - i + 1), nrow = 1)
    V <- a %*% vcov[i:L, i:L] %*% t(a)
    se[i] <- sqrt(V)
  }
  return(se)
}

ss_rev_result <- ss_serevcumsum(test_vcov)
diff_serevcumsum <- max(abs(our_serevcumsum - ss_rev_result))
cat("serevcumsum max diff:", format(diff_serevcumsum, scientific = TRUE), "\n")

if (diff_secumsum == 0 && diff_serevcumsum == 0) {
  cat("PASS: secumsum and serevcumsum are BIT-IDENTICAL to S&S functions\n")
} else {
  cat("Diffs found (should be 0):", diff_secumsum, diff_serevcumsum, "\n")
}

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

results <- data.frame(
  Test = c(
    "secumsum = S&S secumsum",
    "serevcumsum = S&S serevcumsum",
    "DLM = ES on generate_data (coef)",
    "DLM = ES on generate_data (SE)",
    "DLM gamma = plm gamma on BF2017",
    "DLM = ES on BF2017 (coef)",
    "DLM = ES on BF2017 (SE)"
  ),
  Max_Diff = c(
    diff_secumsum,
    diff_serevcumsum,
    max_coef_diff,
    max_se_diff,
    gamma_diff,
    max(coef_diffs),
    max(se_diffs)
  ),
  stringsAsFactors = FALSE
)
results$Pass <- results$Max_Diff < 1e-10
print(results)
cat("\n")
