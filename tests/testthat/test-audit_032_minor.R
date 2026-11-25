# Minor findings of the final pre-release audit of 0.3.2 (F3-F7) and the
# criterion-prefix compatibility note (section 9)

test_that("seeds have their own validator (F3)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  expect_identical(check_seed(20261118), 20261118L)
  expect_identical(check_seed(-5), -5L)
  expect_error(check_seed(2^40), "whole number between")
  expect_error(check_seed(1.5), "whole number between")
  skip_if_not_installed("bootUR")
  expect_s3_class(boot_ur(y, B = 99, seed = 20261118), "ur_test")
})

test_that("an absurd frequency inherited from a ts fails cleanly (F4)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  z <- stats::ts(y, frequency = .Machine$integer.max)
  expect_error(recursive_adf(z), "frequency of the series .* too large")
  expect_error(residual_check(adf(z, lags = 0)), "frequency of the series .* too large")
  # the test itself still works; the residual check is reported as not done
  r <- adf(z, lags = 0)
  expect_s3_class(r, "ur_test")
  expect_match(r$warnings, "Residual check not done", all = FALSE)
  expect_s3_class(recursive_adf(z, frequency = 12), "ur_recursive")
})

test_that("only a full-sample lag choice is reported as look-ahead (F5)", {
  set.seed(1)
  y <- cumsum(stats::rnorm(100))
  look <- function(r) any(grepl("look-ahead", capture.output(print(r))))
  expect_false(look(recursive_adf(y, lags = 0)))
  expect_true(look(recursive_adf(y)))
  expect_false(look(recursive_adf(y, reselect = TRUE)))
  expect_identical(attr(recursive_adf(y, lags = 0), "settings")$lag_mode, "fixed")
})

test_that("the simulation plot accepts graphical arguments (F6)", {
  s <- simulate_df(T = 20, R = 100, seed = 1)
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f); on.exit(unlink(f))
  expect_silent(plot(s, main = "Custom", xlab = "tau", breaks = 20))
  grDevices::dev.off()
})

test_that("a zero bandwidth gives the zero-bandwidth limit gamma0 (F7)", {
  y <- c(1, 0, -1, 0)
  r <- kpss(y, kernel = "qs", bandwidth = "andrews")
  expect_equal(r$statistic, 0.25)
  expect_equal(r$other_statistics[["long_run_variance"]], 0.5)
  for (k in c("bartlett", "parzen")) {
    expect_equal(kpss(y, kernel = k, bandwidth = "andrews")$statistic, 0.25)
  }
  expect_equal(long_run_variance(c(1, -1, 1, -1), "qs", 0)$lrv, 1)
})

test_that("unambiguous criterion abbreviations work again (section 9)", {
  set.seed(2)
  y <- cumsum(stats::rnorm(150))
  expect_identical(adf(y, criterion = "AI")$lag_selection$method, "AIC")
  expect_identical(adf(y, criterion = "ai")$lag_selection$method, "AIC")
  expect_identical(adf(y, criterion = "MA")$lag_selection$method, "MAIC")
  expect_error(adf(y, criterion = "M"), "unambiguous abbreviation")
  expect_error(adf(y, criterion = ""), "must be one of")
})
