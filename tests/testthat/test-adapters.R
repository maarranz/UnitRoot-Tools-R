skip_if_not_installed("urca")

sim_rw <- function(n = 150, seed = 8) {
  set.seed(seed)
  cumsum(stats::rnorm(n))
}

test_that("pp() returns urca's Z_tau with MacKinnon critical values", {
  y <- sim_rw()
  for (trend in c("c", "ct")) for (bw in list("short", "long", 3L)) {
    r <- pp(y, trend = trend, bandwidth = bw)
    u <- urca::ur.pp(y, type = "Z-tau",
                     model = if (trend == "c") "constant" else "trend",
                     lags = if (is.character(bw)) bw else "short",
                     use.lag = if (is.numeric(bw)) bw else NULL)
    expect_s3_class(r, "ur_test")
    expect_equal(r$statistic, unname(u@teststat), tolerance = 1e-12)
    expect_identical(r$lags, u@lag)
    expect_identical(r$nobs, length(y) - 1L)
    expect_equal(r$critical_values, tau_critical_values(trend, r$nobs))
    expect_match(r$engine, "^urca ")
  }
  expect_match(pp(y)$bandwidth_rule, "short")
  expect_identical(pp(y, bandwidth = 2)$bandwidth_rule, "fixed by the user")
})

test_that("kpss() returns urca's statistic with upper-tail decisions", {
  y <- sim_rw()
  for (trend in c("c", "ct")) {
    r <- kpss(y, trend = trend)
    u <- urca::ur.kpss(y, type = if (trend == "c") "mu" else "tau", lags = "short")
    expect_equal(r$statistic, unname(u@teststat), tolerance = 1e-12)
    expect_identical(r$tail, "upper")
    expect_identical(r$alternative, "unit root")
    expect_equal(unname(r$critical_values[c("0.1", "0.05", "0.025", "0.01")]),
                 unname(u@cval[1, ]))
    expect_identical(r$reject, r$statistic > r$critical_values[["0.05"]])
  }
  expect_identical(kpss(y)$statistic_name, "eta_mu")
})

test_that("dfgls() returns urca's DF-GLS statistic and requires lags", {
  y <- sim_rw()
  for (trend in c("c", "ct")) for (k in c(0, 3)) {
    r <- dfgls(y, trend = trend, lags = k)
    u <- urca::ur.ers(y, type = "DF-GLS",
                      model = if (trend == "c") "constant" else "trend", lag.max = k)
    expect_equal(r$statistic, unname(u@teststat), tolerance = 1e-12)
    expect_identical(r$nobs, length(y) - as.integer(k) - 1L)
    expect_identical(r$nobs, length(stats::residuals(u@testreg)))
  }
  expect_error(dfgls(y), "Give the number of lagged differences")
  expect_error(dfgls(y, lags = -1), "non-negative whole number")
})

test_that("adapters check inputs before calling urca", {
  y <- sim_rw()
  expect_error(pp(c(y, NA)), "missing or non-finite")
  expect_error(dfgls(c(NA, y), lags = 1), "missing or non-finite")    # urca would drop it
  expect_error(pp(y, trend = "n"), "not available here")
  expect_error(kpss(y, bandwidth = "andrews"), "not available in this version")
  expect_error(kpss(y, bandwidth = 1.5), "non-negative whole number")
})

test_that("adapter results print like every other test", {
  out <- capture.output(print(kpss(sim_rw(), trend = "ct")))
  expect_match(out, "KPSS test of stationarity \\(eta_tau\\)", all = FALSE)
  expect_match(out, "H0: stationary around a linear trend", all = FALSE)
  expect_match(out, "Bandwidth: \\d+ \\(short", all = FALSE)
})
