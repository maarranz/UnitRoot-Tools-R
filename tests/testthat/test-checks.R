# Robustness checks added after the Phase 1 audit (findings M1, M2, M3, m1-m3)

set.seed(31)
y <- cumsum(stats::rnorm(150))

test_that("M1: arrays and matrix columns are not flattened into a series", {
  expect_error(check_series(array(stats::rnorm(60), c(10, 3, 2))), "3 dimensions")
  df <- data.frame(date = as.Date("2020-01-01") + 0:9)
  df$m <- matrix(stats::rnorm(20), 10)
  expect_error(check_series(df, var = "m"), "matrix or array")
  df$v <- stats::rnorm(10)
  expect_identical(check_series(df), df$v)            # the matrix column is not a candidate
  expect_error(adf(array(stats::rnorm(60), c(10, 3, 2))), "3 dimensions")
})

test_that("M2: orders are never truncated or coerced", {
  for (bad in list(2.9, -0.5, -1, "3", NA, Inf, c(1, 2), TRUE)) {
    expect_error(select_lags(y, max_lags = bad), "`max_lags` must be a single non-negative whole number")
    expect_error(adf(y, lags = bad), "`lags` must be a single non-negative whole number")
    expect_error(ur_strategy(y, lags = bad), "`lags` must be a single non-negative whole number")
  }
  expect_error(dfgls(y, lags = 1.5), "`lags` must be a single")
  expect_error(pp(y, bandwidth = 2.5), "`bandwidth` must be a single")
  expect_identical(select_lags(y, max_lags = 3)$kmax, 3L)
  expect_identical(select_lags(y, max_lags = 3L)$kmax, 3L)
})

test_that("m1: every regression keeps two residual degrees of freedom", {
  expect_error(qr_fit(cbind(1, 1:3), c(1, 3, 2)), "two residual degrees of freedom")
  expect_identical(qr_fit(cbind(1, 1:4), c(1, 3, 2, 5))$df, 2L)
})

test_that("m2: allow_no_constant must be TRUE or FALSE", {
  expect_error(ur_strategy(y, "conservative", allow_no_constant = 2, lags = 1), "TRUE or FALSE")
  expect_error(ur_strategy(y, "conservative", allow_no_constant = NA, lags = 1), "TRUE or FALSE")
  expect_error(ur_strategy(y, "conservative", allow_no_constant = c(TRUE, FALSE), lags = 1),
               "TRUE or FALSE")
})

test_that("m3: levels are checked up front, with the available ones listed", {
  expect_error(adf(y, level = 0.025), "one of 0.01, 0.05, 0.1")
  expect_error(pp(y, level = 0.2), "one of 0.01, 0.05, 0.1")
  expect_error(dfgls(y, lags = 1, level = "0.05"), "`level` must be one of")
  expect_error(ur_strategy(y, level = 0.025), "one of 0.01, 0.05, 0.1")
  expect_true(kpss(y, level = 0.025)$level == 0.025)  # KPSS tabulates 2.5%
  expect_error(kpss(y, level = 0.2), "0.01, 0.025, 0.05, 0.1")
})
