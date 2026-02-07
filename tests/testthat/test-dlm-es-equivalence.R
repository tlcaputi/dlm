# test-dlm-es-equivalence.R
# Tests that DLM betas and SEs match event-study betas and SEs exactly.

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

  # SE(g1+g2) = sqrt([1,1] %*% V3[1:2,1:2] %*% [1,1]') = sqrt(1 + 1 + 2*0.5) = sqrt(3)
  # SE(g1+g2+g3) = sqrt([1,1,1] %*% V3 %*% [1,1,1]') = sqrt(1+1+1 + 2*(0.5+0+0.5)) = sqrt(5)
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
  # Rev cumsum:
  # i=1: sum from 1 to 3 -> same as secumsum full = sqrt(5)
  # i=2: sum from 2 to 3 -> sqrt([1,1] %*% V3[2:3,2:3] %*% [1,1]') = sqrt(1+1+2*0.5) = sqrt(3)
  # i=3: just g3 -> sqrt(1) = 1
  result3 <- serevcumsum(V3)
  expect_equal(result3, c(sqrt(5), sqrt(3), 1))
})

test_that("DLM betas match event-study betas exactly", {
  skip_if_not_installed("fixest")

  df <- generate_data(seed = 42)

  from_rt <- -5
  to_rt <- 5
  ref_period <- -1

  # Run DLM
  dlm_results <- distributed_lags_models(
    data          = df,
    exposure_data = df,
    from_rt       = from_rt,
    to_rt         = to_rt,
    outcomes      = "outcome",
    exposure      = "post",
    unit          = "group",
    time          = "time",
    ref_period    = ref_period
  )

  dlm_betas <- dlm_results[[1]]$betas

  # Run equivalent event study with fixest::i() and bin
  es_model <- fixest::feols(
    outcome ~ i(years_to_treatment, treat,
                ref = c(ref_period, -1000),
                bin = list(
                  "-5+" = ~ x <= -5,
                  "5+"  = ~ x >= 5
                )) | group + time,
    data = df,
    cluster = ~group
  )

  es_ct <- as.data.frame(es_model$coeftable)
  es_betas <- es_ct[, 1]  # coefficients
  es_ses   <- es_ct[, 2]  # standard errors

  # DLM betas should match ES betas
  expect_equal(
    dlm_betas$coef,
    es_betas,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )

  # DLM SEs should match ES SEs
  expect_equal(
    dlm_betas$se,
    es_ses,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("DLM recovers true effect of -3", {
  df <- generate_data(seed = 123, n_groups = 26^2, n_times = 50, treat_prob = 0.5)

  dlm_results <- distributed_lags_models(
    data          = df,
    exposure_data = df,
    from_rt       = -5,
    to_rt         = 5,
    outcomes      = "outcome",
    exposure      = "post",
    unit          = "group",
    time          = "time",
    ref_period    = -1
  )

  post_betas <- dlm_results[[1]]$betas
  post_coefs <- post_betas$coef[post_betas$time_to_event >= 0]

  # All post-treatment coefficients should be close to -3
  for (b in post_coefs) {
    expect_true(abs(b - (-3)) < 1,
                info = paste("Post-treatment beta", b, "not within 1 of true effect -3"))
  }
})

test_that("DLM works with single lead", {
  df <- generate_data(seed = 99)

  results <- distributed_lags_models(
    data          = df,
    exposure_data = df,
    from_rt       = -1,
    to_rt         = 3,
    outcomes      = "outcome",
    exposure      = "post",
    unit          = "group",
    time          = "time",
    ref_period    = -1
  )

  expect_false(is.null(results))
  expect_equal(nrow(results[[1]]$betas), 3)  # periods 0, 1, 2, 3 minus ref = 3 betas
})

test_that("DLM works with single lag", {
  df <- generate_data(seed = 99)

  results <- distributed_lags_models(
    data          = df,
    exposure_data = df,
    from_rt       = -3,
    to_rt         = 1,
    outcomes      = "outcome",
    exposure      = "post",
    unit          = "group",
    time          = "time",
    ref_period    = -1
  )

  expect_false(is.null(results))
  # Periods: -3, -2, 0, 1 (ref=-1 excluded) = 4 betas
  expect_equal(nrow(results[[1]]$betas), 4)
})

test_that("DLM works with from_rt == ref_period", {
  df <- generate_data(seed = 77)

  results <- distributed_lags_models(
    data          = df,
    exposure_data = df,
    from_rt       = -1,
    to_rt         = 5,
    outcomes      = "outcome",
    exposure      = "post",
    unit          = "group",
    time          = "time",
    ref_period    = -1
  )

  expect_false(is.null(results))
  # Periods: 0, 1, 2, 3, 4, 5 = 6 betas
  expect_equal(nrow(results[[1]]$betas), 6)
})
