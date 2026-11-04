skip_if_not_installed("urca")

sim_rw <- function(n = 150, seed = 8) {
  set.seed(seed)
  cumsum(stats::rnorm(n))
}

test_that("pp(): internal Z_tau, bandwidth as urca, MacKinnon calibration (D12)", {
  y <- sim_rw()
  for (trend in c("c", "ct")) for (bw in list("short", "long", 3L)) {
    r <- pp(y, trend = trend, bandwidth = bw)
    u <- urca::ur.pp(y, type = "Z-tau",
                     model = if (trend == "c") "constant" else "trend",
                     lags = if (is.character(bw)) bw else "short",
                     use.lag = if (is.numeric(bw)) bw else NULL)
    expect_s3_class(r, "ur_test")
    expect_identical(r$lags, u@lag)                      # same bandwidth rule
    expect_equal(r$statistic, unname(u@teststat), tolerance = 0.05)  # O(1/T) apart
    expect_identical(r$nobs, length(y) - 1L)
    expect_identical(r$sample, c(2L, length(y)))
    expect_equal(r$critical_values, tau_critical_values(trend, r$nobs))
    expect_equal(r$p_value, tau_p_value(r$statistic, trend, r$nobs))
    expect_match(r$engine, "^internal ")
  }
  expect_match(pp(y)$bandwidth_rule, "short rule, trunc\\(4 \\(\\(T-1\\)/100\\)", fixed = FALSE)
  expect_identical(pp(y, bandwidth = 2)$bandwidth_rule,
                   "Bartlett kernel, 2 autocovariances (fixed by the user)")
  expect_error(pp(y[1:10], bandwidth = 9), "too large")
  # with bandwidth 0 there is no correction: Z_tau is the Dickey-Fuller t ratio
  expect_equal(pp(y, bandwidth = 0)$statistic,
               adf(y, trend = "c", lags = 0)$statistic, tolerance = 1e-12)
})

test_that("kpss(): internal statistic equals urca's; exact asymptotic critical values (D14)", {
  y <- sim_rw()
  for (trend in c("c", "ct")) {
    r <- kpss(y, trend = trend)
    u <- urca::ur.kpss(y, type = if (trend == "c") "mu" else "tau", lags = "short")
    expect_equal(r$statistic, unname(u@teststat), tolerance = 1e-12)
    expect_identical(r$tail, "upper")
    expect_identical(r$alternative, "unit root")
    # Kwiatkowski et al.'s simulated table (used by urca) is close but not exact
    expect_equal(unname(r$critical_values[c("0.1", "0.05", "0.025", "0.01")]),
                 unname(u@cval[1, ]), tolerance = 0.015)
    expect_match(r$engine, "^internal ")
    expect_true(r$p_value > 0 && r$p_value < 1)
    expect_identical(r$reject, r$statistic > r$critical_values[["0.05"]])
  }
  expect_identical(kpss(y)$statistic_name, "eta_mu")
})

test_that("dfgls() is computed internally and agrees with urca's DF-GLS statistic", {
  y <- sim_rw()
  expect_match(dfgls(y, lags = 1)$engine, "^internal ")
  for (trend in c("c", "ct")) for (k in c(0, 3, 8)) {
    r <- dfgls(y, trend = trend, lags = k)
    u <- urca::ur.ers(y, type = "DF-GLS",
                      model = if (trend == "c") "constant" else "trend", lag.max = k)
    expect_equal(r$statistic, unname(u@teststat), tolerance = 1e-12)
    expect_identical(r$nobs, length(y) - as.integer(k) - 1L)
    expect_identical(r$nobs, length(stats::residuals(u@testreg)))
  }
  expect_error(dfgls(y, lags = -1), "non-negative whole number")
  # audit N1: two residual degrees of freedom, as for the native regressions
  expect_error(dfgls(c(1, 3, 2), trend = "c", lags = 0), "two residual degrees of freedom")
  expect_error(dfgls(cumsum(rnorm(9)), trend = "c", lags = 3), "two residual degrees of freedom")
  expect_s3_class(dfgls(cumsum(rnorm(10)), trend = "c", lags = 3), "ur_test")
  expect_error(ng_perron(cumsum(rnorm(9)), trend = "c", lags = 3), "two residual degrees of freedom")
  expect_error(dfgls(y, lags = 2, max_lags = 4), "either `lags`")
})

test_that("dfgls() selects lags by Perron-Qu: MAIC on the OLS-detrended series (D9)", {
  set.seed(5)
  y <- cumsum(stats::arima.sim(list(ar = 0.6), n = 200))
  for (trend in c("c", "ct")) for (crit in c("MAIC", "MBIC")) {
    r <- dfgls(y, trend = trend, criterion = crit)
    sel <- select_lags(y, trend = trend, criterion = crit)
    expect_identical(r$lags, sel$lags)
    expect_equal(r$lag_selection$table, sel$table)
    expect_identical(r$lag_selection$kmax, sel$kmax)
    expect_match(r$lag_selection$detrending, "Perron and Qu")
    expect_equal(r$statistic, dfgls(y, trend = trend, lags = r$lags)$statistic)
    expect_identical(r$nobs, length(y) - r$lags - 1L)     # maximal sample
  }
  expect_identical(dfgls(y, max_lags = 3)$lag_selection$kmax, 3L)
  expect_match(capture.output(print(dfgls(y))), "Lags: \\d+ \\(MAIC, kmax = 14\\)", all = FALSE)
})

test_that("inputs are checked before any computation", {
  y <- sim_rw()
  expect_error(pp(c(y, NA)), "missing or non-finite")
  expect_error(dfgls(c(NA, y), lags = 1), "missing or non-finite")    # urca would drop it
  expect_error(pp(y, trend = "n"), "not available here")
  expect_error(kpss(y, bandwidth = "automatic"), "must be \"short\"")
  expect_error(kpss(y, kernel = "tukey"), "`kernel` must be one of")
  expect_error(kpss(y, bandwidth = 1.5), "non-negative whole number")
})

test_that("adapter results print like every other test", {
  out <- capture.output(print(kpss(sim_rw(), trend = "ct")))
  expect_match(out, "KPSS test of stationarity \\(eta_tau\\)", all = FALSE)
  expect_match(out, "H0: stationary around a linear trend", all = FALSE)
  expect_match(out, "Long-run variance: Bartlett kernel, \\d+ autocovariances \\(short rule", all = FALSE)
})
