sim_series <- function(n = 150, seed = 21) {
  set.seed(seed)
  cumsum(0.2 + stats::arima.sim(list(ar = 0.4), n = n))
}

test_that("fixed lags reproduce urca::ur.df (tau and phi statistics)", {
  skip_if_not_installed("urca")
  y <- sim_series()
  urca_type <- c(n = "none", c = "drift", ct = "trend")
  for (trend in c("n", "c", "ct")) for (k in c(0, 1, 4)) {
    r <- adf(y, trend = trend, lags = k)
    u <- urca::ur.df(y, type = urca_type[[trend]], lags = k)
    stats_u <- u@teststat[1, ]
    expect_equal(r$statistic, unname(stats_u[1]), tolerance = 1e-8)
    if (trend != "n") {
      expect_equal(unname(r$other_statistics), unname(stats_u[-1]), tolerance = 1e-8)
    }
    expect_identical(r$nobs, as.integer(length(y) - k - 1))   # maximal sample, as in urca
  }
})

test_that("tau and its critical values follow the specification", {
  y <- sim_series()
  r <- adf(y, trend = "ct", lags = 2)
  expect_s3_class(r, "ur_test")
  expect_identical(r$statistic_name, "tau3")
  expect_identical(r$tail, "lower")
  expect_equal(r$critical_values, tau_critical_values("ct", r$nobs))
  expect_identical(r$calibration$source, "MacKinnon (2010) response surface")
  expect_true(is.na(r$p_value))
  expect_named(r$other_statistics, c("phi2", "phi3"))
  expect_named(r$other_critical_values, c("phi2", "phi3"))
  expect_identical(adf(y, trend = "n", lags = 0)$calibration$source,
                   "MacKinnon (1996) response surface")
  expect_identical(r$reject, r$statistic < r$critical_values[["0.05"]])
})

test_that("lag selection: common sample by default, maximal on request", {
  y <- sim_series(n = 200)
  sel <- select_lags(y, trend = "c")
  r1 <- adf(y, trend = "c")
  expect_identical(r1$lags, sel$lags)
  expect_identical(r1$nobs, sel$nobs)
  expect_identical(r1$lag_selection$final_sample, "common")
  expect_equal(r1$lag_selection$table, sel$table)

  r2 <- adf(y, trend = "c", final_sample = "maximal")
  expect_identical(r2$nobs, length(y) - r2$lags - 1L)
  expect_equal(r2$statistic, adf(y, trend = "c", lags = r2$lags)$statistic)
})

test_that("argument checks", {
  y <- sim_series()
  expect_error(adf(y, lags = 2, max_lags = 4), "either `lags`")
  expect_error(adf(y, lags = 1.5), "non-negative integer")
  expect_error(adf(y, trend = "ctt"), "Unknown deterministic")
  expect_error(adf(c(y, NA)), "missing or non-finite")
})

test_that("small samples are flagged", {
  set.seed(2)
  r <- adf(cumsum(rnorm(22)), trend = "c", lags = 0)
  expect_match(r$warnings, "below 25", all = FALSE)
  expect_match(r$warnings, "extrapolation", all = FALSE)
})

test_that("works with ts and tidy inputs", {
  y <- sim_series()
  r_num <- adf(y, trend = "c", lags = 1)
  expect_equal(adf(ts(y, start = 1990, frequency = 4), trend = "c", lags = 1)$statistic,
               r_num$statistic)
  df <- data.frame(date = as.Date("2000-01-01") + seq_along(y), gdp = y)
  expect_equal(adf(df, trend = "c", lags = 1)$statistic, r_num$statistic)
})

test_that("print and summary show the test, phi statistics and lag table", {
  y <- sim_series()
  r <- adf(y, trend = "ct")
  out <- capture.output(print(r))
  expect_match(out, "Augmented Dickey-Fuller test \\(tau3\\)", all = FALSE)
  expect_match(out, "Lags: \\d+ \\(MAIC, kmax = \\d+\\)", all = FALSE)
  s <- capture.output(print(summary(r)))
  expect_match(s, "Accompanying statistics", all = FALSE)
  expect_match(s, "phi3", all = FALSE)
  expect_match(s, "Dickey and Fuller \\(1981\\), table row T = 100", all = FALSE)
  expect_match(s, "test regression on the common sample", all = FALSE)
  expect_match(s, "y_lag1", all = FALSE)
})
