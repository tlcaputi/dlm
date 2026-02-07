# test-dlm-es-equivalence.R
# Tests that DLM betas and SEs match event-study betas and SEs exactly,
# confirming the Schmidheiny & Siegloch (2023, JAE) equivalence result.

test_that("secumsum matches known analytic values", {
  # For a 2x2 identity matrix, cumulative sums should give:
  # SE(gamma1) = sqrt(1) = 1
  # SE(gamma1 + gamma2) = sqrt(1 + 1 + 2*0) = sqrt(2)
  V <- diag(2)
  result <- secumsum(V)
  expect_equal(result, c(1, sqrt(2)))

  # For a single element
  expect_equal(secumsum(matrix(4)), 2)

  # For a 3x3 with off-diagonals
  V3 <- matrix(c(1, 0.5, 0,
                  0.5, 1, 0.5,
                  0, 0.5, 1), nrow = 3)
  # SE(g1) = sqrt(1) = 1
  # SE(g1+g2) = sqrt([1,1] V[1:2,1:2] [1,1]') = sqrt(1+1+2*0.5) = sqrt(3)
  # SE(g1+g2+g3) = sqrt([1,1,1] V [1,1,1]') = sqrt(1+1+1+2*(0.5+0+0.5)) = sqrt(5)
  result3 <- secumsum(V3)
  expect_equal(result3, c(1, sqrt(3), sqrt(5)))
})

test_that("serevcumsum matches known analytic values", {
  V <- diag(2)
  result <- serevcumsum(V)
  # Reverse cumulative sum: SE(g1+g2) at position 1, SE(g2) at position 2
  expect_equal(result, c(sqrt(2), 1))

  expect_equal(serevcumsum(matrix(4)), 2)

  V3 <- matrix(c(1, 0.5, 0,
                  0.5, 1, 0.5,
                  0, 0.5, 1), nrow = 3)
  # i=1: sum from 1 to 3 -> sqrt(5)
  # i=2: sum from 2 to 3 -> sqrt(3)
  # i=3: just g3 -> sqrt(1) = 1
  result3 <- serevcumsum(V3)
  expect_equal(result3, c(sqrt(5), sqrt(3), 1))
})

# Helper: compute the DLM effective time window so we can restrict the ES
# to the same sample. The DLM creates lead(abs(from_rt)-1) through lag(to_rt),
# so it drops observations where any lead/lag is NA:
#   - lag(x, k) is NA for the first k periods -> need time > to_rt
#   - lead(x, k) is NA for the last k periods -> need time <= max_time - (abs(from_rt)-1)
dlm_time_window <- function(n_times, from_rt, to_rt) {
  min_t <- to_rt + 1
  max_t <- n_times - (abs(from_rt) - 1)
  c(min_t, max_t)
}

test_that("DLM betas and SEs match ES betas and SEs exactly (from_rt=-3, to_rt=3)", {
  skip_if_not_installed("fixest")

  df <- generate_data(seed = 42)
  from_rt <- -3; to_rt <- 3; ref_period <- -1

  # Run DLM
  dlm_results <- distributed_lags_models(
    data = df, exposure_data = df,
    from_rt = from_rt, to_rt = to_rt,
    outcomes = "outcome", exposure = "post",
    unit = "group", time = "time",
    ref_period = ref_period
  )
  dlm_betas <- dlm_results[[1]]$betas
  dlm_N <- stats::nobs(dlm_results[[1]]$model)

  # Restrict ES to same time window (critical for equivalence)
  tw <- dlm_time_window(20, from_rt, to_rt)
  df_r <- df[df$time >= tw[1] & df$time <= tw[2], ]
  expect_equal(nrow(df_r), dlm_N)

  # Run equivalent ES with fixest::i() and bin
  es_model <- fixest::feols(
    outcome ~ i(years_to_treatment, treat,
                ref = c(ref_period, -1000),
                bin = list(
                  "-3+" = ~ x <= -3,
                  "3+"  = ~ x >= 3
                )) | group + time,
    data = df_r, cluster = ~group
  )
  es_ct <- as.data.frame(es_model$coeftable)

  # Coefficients must match to machine precision
  expect_equal(dlm_betas$coef, es_ct[, 1], tolerance = 1e-10, ignore_attr = TRUE)

  # SEs must match to machine precision
  expect_equal(dlm_betas$se, es_ct[, 2], tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("DLM betas and SEs match ES exactly (from_rt=-5, to_rt=5, wide treatment)", {
  skip_if_not_installed("fixest")

  # Use wider treatment timing to avoid collinearity at long leads
  set.seed(42)
  n_groups <- 200; n_times <- 30
  groups <- paste0("g", 1:n_groups)
  panel <- expand.grid(group = groups, time = 1:n_times)
  treatments <- data.frame(group = groups) |>
    dplyr::sample_frac(0.5) |>
    dplyr::mutate(treatment_time = sample(5:15, size = dplyr::n(), replace = TRUE))
  df <- merge(panel, treatments, by = "group", all.x = TRUE) |>
    dplyr::mutate(
      treat = as.numeric(!is.na(treatment_time)),
      years_to_treatment = ifelse(treat == 0, -1000, time - treatment_time),
      post = ifelse(treat == 0, 0, as.numeric(time >= treatment_time)),
      outcome = stats::rnorm(dplyr::n(), sd = 5) + ifelse(years_to_treatment >= 0, -3, 0)
    ) |> dplyr::arrange(group, time)

  from_rt <- -5; to_rt <- 5; ref_period <- -1

  dlm_results <- distributed_lags_models(
    data = df, exposure_data = df,
    from_rt = from_rt, to_rt = to_rt,
    outcomes = "outcome", exposure = "post",
    unit = "group", time = "time",
    ref_period = ref_period
  )
  dlm_betas <- dlm_results[[1]]$betas
  dlm_N <- stats::nobs(dlm_results[[1]]$model)

  tw <- dlm_time_window(n_times, from_rt, to_rt)
  df_r <- df[df$time >= tw[1] & df$time <= tw[2], ]
  expect_equal(nrow(df_r), dlm_N)

  es_model <- fixest::feols(
    outcome ~ i(years_to_treatment, treat,
                ref = c(ref_period, -1000),
                bin = list("-5+" = ~ x <= -5, "5+" = ~ x >= 5)) | group + time,
    data = df_r, cluster = ~group
  )
  es_ct <- as.data.frame(es_model$coeftable)

  expect_equal(dlm_betas$coef, es_ct[, 1], tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(dlm_betas$se, es_ct[, 2], tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("DLM recovers true effect of -3", {
  # Use many groups for a precise estimate
  df <- generate_data(seed = 123, n_groups = 2000, n_times = 40, treat_prob = 0.5)

  dlm_results <- distributed_lags_models(
    data = df, exposure_data = df,
    from_rt = -3, to_rt = 3,
    outcomes = "outcome", exposure = "post",
    unit = "group", time = "time",
    ref_period = -1
  )

  post_betas <- dlm_results[[1]]$betas
  post_coefs <- post_betas$coef[post_betas$time_to_event >= 0]

  # All post-treatment coefficients should be within 1 of -3
  for (b in post_coefs) {
    expect_true(abs(b - (-3)) < 1,
                info = paste("Post-treatment beta", round(b, 3), "not within 1 of true effect -3"))
  }
})

test_that("DLM works with from_rt == ref_period", {
  df <- generate_data(seed = 77)

  results <- distributed_lags_models(
    data = df, exposure_data = df,
    from_rt = -1, to_rt = 3,
    outcomes = "outcome", exposure = "post",
    unit = "group", time = "time",
    ref_period = -1
  )

  expect_false(is.null(results))
  # Periods: 0, 1, 2, 3 = 4 betas
  expect_equal(nrow(results[[1]]$betas), 4)
})

test_that("generate_data does not leak RNG state", {
  set.seed(999)
  x_before <- stats::rnorm(5)

  set.seed(999)
  df <- generate_data(seed = 42)
  x_after <- stats::rnorm(5)

  expect_equal(x_before, x_after)
})
