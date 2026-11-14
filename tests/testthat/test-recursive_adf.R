# Recursive and rolling ADF statistics (spec: Decision D20)

rw_m <- function(n = 240, freq = 12) {
  set.seed(7)
  stats::ts(cumsum(stats::rnorm(n)), start = 2000, frequency = freq)
}

test_that("each window is the ADF test on that subsample, with fixed lags", {
  y <- rw_m()
  rr <- recursive_adf(y, trend = "c")
  k <- adf(y, trend = "c")$lags
  expect_true(all(rr$lags == k))
  last <- adf(y, trend = "c", lags = k)
  expect_equal(rr$statistic[nrow(rr)], last$statistic)
  expect_equal(rr$critical_value[nrow(rr)], last$critical_values[["0.05"]])
  expect_equal(rr$p_value[nrow(rr)], last$p_value)
  # a window in the middle
  i <- 50
  mid <- adf(as.numeric(y)[1:(60 + i - 1)], trend = "c", lags = k)
  expect_equal(rr$statistic[i], mid$statistic)
  expect_identical(rr$reject, rr$statistic < rr$critical_value)
})

test_that("default sizes follow the frequency, with floors", {
  y <- rw_m()
  rr <- recursive_adf(y, lags = 1)
  expect_identical(nrow(rr), 240L - 60L + 1L)     # monthly: first sample 5s = 60
  ro <- recursive_adf(y, lags = 1, window = "rolling")
  expect_identical(nrow(ro), 240L - 120L + 1L)    # monthly: width 10s = 120
  expect_true(all(ro$nobs == 120L - 2L))
  q <- rw_m(freq = 4)
  expect_identical(attr(recursive_adf(q, lags = 1), "settings")$size, 30L)   # max(20, 30)
  expect_identical(attr(recursive_adf(q, lags = 1, window = "rolling"), "settings")$size, 40L)
  v <- as.numeric(y)
  rv <- recursive_adf(v, lags = 1)
  expect_identical(attr(rv, "settings")$size, 30L)
  expect_match(attr(rv, "settings")$notes, "frequency unknown", all = FALSE)
  expect_identical(attr(recursive_adf(v, lags = 1, frequency = 12), "settings")$size, 60L)
  expect_identical(nrow(recursive_adf(y, lags = 1, min_size = 100)), 141L)
})

test_that("time points come from the series", {
  y <- rw_m()
  rr <- recursive_adf(y, lags = 1)
  expect_equal(rr$end[1], 2000 + 59 / 12)
  expect_equal(rr$start[1], 2000)
  d <- data.frame(date = seq(as.Date("2000-01-01"), by = "month", length.out = 240), x = as.numeric(y))
  rd <- recursive_adf(d, lags = 1, window = "rolling")
  expect_s3_class(rd$end, "Date")
  expect_identical(rd$start[2], as.Date("2000-02-01"))
})

test_that("lags can be re-selected in each window", {
  y <- rw_m()
  rs <- recursive_adf(y, reselect = TRUE)
  expect_match(attr(rs, "settings")$lag_choice, "selected in each window")
  j <- 30
  sub <- as.numeric(y)[1:(60 + j - 1)]
  expect_identical(rs$lags[j], adf(sub, trend = "c")$lags)
  expect_error(recursive_adf(y, lags = 2, reselect = TRUE), "needs lag selection")
})

test_that("invalid sizes stop with a message", {
  y <- rw_m(60)
  expect_error(recursive_adf(y, lags = 1, window = "rolling"), "longer than the series")
  expect_error(recursive_adf(y, lags = 12, min_size = 25, trend = "ct"), "too small")
})

test_that("print and plot work", {
  rr <- recursive_adf(rw_m(), lags = 1)
  out <- capture.output(print(rr))
  expect_match(out, "Recursive ADF test", all = FALSE)
  expect_match(out, "pointwise", all = FALSE)
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f); on.exit(unlink(f))
  expect_invisible(plot(rr))
  grDevices::dev.off()
})
