# Regression tests for the audit of 0.3.0 (14 November 2026)

test_that("recursive_adf() respects max_lags in every window when reselecting (H1)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  expect_identical(unique(recursive_adf(y, max_lags = 0, reselect = TRUE)$lags), 0L)
  expect_true(all(recursive_adf(y, max_lags = 2, reselect = TRUE)$lags <= 2L))
  # a cap larger than a small window allows is reduced to the window's limit
  r <- recursive_adf(y, max_lags = 12, reselect = TRUE)
  expect_true(all(r$lags <= 12L))
  expect_match(attr(r, "settings")$lag_choice, "kmax = 12")
})

test_that("KPSS is invariant to the units of the series (H2)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  d <- diff(y)
  base <- kpss(d)$statistic
  for (s in c(1e-300, 1e-12, 1e12, 1e300, -3)) {
    expect_equal(kpss(s * d)$statistic, base, tolerance = 1e-10)
    expect_equal(kpss(s * y, "ct", kernel = "qs", bandwidth = "andrews")$statistic,
                 kpss(y, "ct", kernel = "qs", bandwidth = "andrews")$statistic, tolerance = 1e-10)
  }
  # degenerate series are still refused, whatever the units
  for (s in c(1e-12, 1, 1e12)) {
    expect_error(kpss(rep(s, 100)), "the series is constant")
    expect_error(kpss(s * (3 + 0.5 * seq_len(100)), trend = "ct"), "exact linear trend")
  }
  expect_error(kpss(rep(0, 100)), "the series is constant")
})

test_that("recursive_adf() states the look-ahead and nominal pointwise references (M1)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  out <- capture.output(recursive_adf(y, lags = 1))
  expect_match(out, "nominal pointwise references", all = FALSE)
  expect_match(out, "look-ahead", all = FALSE)
  expect_false(any(grepl("look-ahead", capture.output(recursive_adf(y, reselect = TRUE)))))
})

test_that("an unavailable residual diagnostic is NA, not white noise (M2)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  r <- residual_check(adf(y, lags = 0), bg_orders = 0)
  expect_identical(nrow(r$breusch_godfrey), 0L)
  expect_identical(r$white_noise, NA)
  expect_match(capture.output(print(r)), "Residual check unavailable", all = FALSE)
  expect_true(isTRUE(residual_check(adf(y, lags = 0))$white_noise) ||
              isFALSE(residual_check(adf(y, lags = 0))$white_noise))
})

test_that("extreme whole numbers are refused before any arithmetic (L1)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  big <- .Machine$integer.max
  expect_error(dfgls(y, lags = big), "too large")
  expect_error(ng_perron(y, lags = big), "too large")
  expect_error(adf(y, lags = big), "too large")
  expect_error(recursive_adf(y, frequency = big), "too large")
  expect_error(residual_check(adf(y, lags = 0), frequency = big), "too large")
  expect_error(adf(y, max_lags = 2^40), "too large")
  # ordinary large values still give the usual messages
  expect_error(adf(y, lags = 60), "Not enough observations|too large")
})

test_that("plot arguments override the defaults (L2)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f); on.exit(unlink(f))
  expect_invisible(plot(recursive_adf(y, lags = 0), main = "Custom", xlab = "Year",
                        ylab = "tau", ylim = c(-5, 1), type = "b", cex.main = 1.2))
  expect_invisible(plot(residual_check(adf(y, lags = 0)), main = "Custom", xlab = "k"))
  grDevices::dev.off()
})
