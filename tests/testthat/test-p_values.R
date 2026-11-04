# p-values from the package's own calibrations (audit C1, 4 Nov 2026): the
# distribution function is interpolated monotonically between the simulated
# quantiles and bounded beyond them, never extrapolated.

calibrations <- function() {
  list(dfgls_c = dfgls_rs$c, dfgls_ct = dfgls_rs$ct,
       phi1 = phi_rs$phi1, phi2 = phi_rs$phi2, phi3 = phi_rs$phi3,
       rho_c = rho_rs$c, rho_ct = rho_rs$ct,
       MZa_c = gls_rs$MZa$c, MZt_ct = gls_rs$MZt$ct, MSB_c = gls_rs$MSB$c,
       MPT_ct = gls_rs$MPT$ct, PT_c = gls_rs$PT$c)
}

test_that("p-values are monotone, bounded, and exact at the simulated quantiles", {
  for (nm in names(calibrations())) {
    rs <- calibrations()[[nm]]
    for (T in c(30, 100, 250, 1000, 5000)) for (k in c(0L, 2L)) {
      q <- unname(rs_quantiles(rs, T, k))
      s <- seq(q[1] - 10 * diff(range(q)), q[length(q)] + 10 * diff(range(q)), length.out = 400)
      for (tl in c("lower", "upper")) {
        p <- vapply(s, function(x) rs_p_value(rs, x, T, k, tail = tl), numeric(1))
        d <- diff(p)
        if (tl == "lower") expect_true(all(d >= -1e-12), label = paste(nm, T, k, tl))
        else expect_true(all(d <= 1e-12), label = paste(nm, T, k, tl))
        expect_true(all(p >= 0.001 - 1e-12 & p <= 0.999 + 1e-12), label = paste(nm, T, k, tl))
      }
      # exact at the 1%, 5% and 10% quantiles
      j <- vapply(c(0.01, 0.05, 0.10), function(a) which.min(abs(rs$probs - a)), 1L)
      pv <- vapply(q[j], function(x) rs_p_value(rs, x, T, k), numeric(1))
      expect_equal(pv, c(0.01, 0.05, 0.10), tolerance = 1e-10, label = paste(nm, T, k))
    }
  }
})

test_that("stationary data get small p-values, consistent with rejection (audit C1)", {
  set.seed(1)
  y <- stats::rnorm(1000)
  a <- dfgls(y, trend = "c", lags = 0)
  b <- ng_perron(y, trend = "c", lags = 0)
  for (r in list(a, b)) {
    expect_true(r$reject)
    expect_equal(r$p_value, 0.001)
    expect_match(r$calibration$note, "beyond the simulated quantiles")
    expect_match(capture.output(print(r)), "p-value: +< 0.001", all = FALSE)
  }
  set.seed(1)
  z <- ng_perron(stats::rnorm(250), "c", lags = 0)
  expect_true(z$reject)
  expect_lt(z$p_value, 0.01)
  # every own calibration: decision and p-value agree at the test's level
  set.seed(5)
  for (n in c(60, 250)) {
    x <- stats::arima.sim(list(ar = 0.5), n = n)
    for (r in list(dfgls(x, "ct"), ng_perron(x, "c"), ng_perron(x, "ct", statistic = "MZa"),
                   ers_pt(x, "c"), ng_perron(x, "c", statistic = "MSB"))) {
      expect_identical(r$reject, r$p_value < r$level, label = paste(r$method, n))
    }
  }
})

test_that("phi and Z_alpha p-values are bounded far in the tail", {
  set.seed(2)
  x <- stats::rnorm(500)
  r <- adf(x, trend = "ct", lags = 0)
  expect_equal(attr(r$other_critical_values$phi3, "p_value"), 0.001)
  p <- pp(x, trend = "c")
  expect_equal(attr(p$other_critical_values$Z_alpha, "p_value"), 0.001)
})
