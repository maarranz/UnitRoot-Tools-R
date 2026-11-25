# Regression tests for the final pre-release audit of 0.3.2 (18 November 2026)

test_that("KPSS variance fields are in the units of the series (major 1)", {
  set.seed(1)
  d <- diff(cumsum(stats::rnorm(200)))
  r1 <- kpss(d)
  r10 <- kpss(10 * d)
  expect_equal(r1$other_statistics[["gamma0"]], mean((d - mean(d))^2), tolerance = 1e-12)
  expect_equal(r10$other_statistics[["gamma0"]], 100 * r1$other_statistics[["gamma0"]], tolerance = 1e-12)
  expect_equal(r10$other_statistics[["long_run_variance"]],
               100 * r1$other_statistics[["long_run_variance"]], tolerance = 1e-12)
  expect_equal(r10$statistic, r1$statistic, tolerance = 1e-12)
  # the long-run variance is the Bartlett estimate in the original units
  l <- r1$lags
  n <- length(d)
  g <- vapply(0:l, function(j) sum((d - mean(d))[(j + 1):n] * (d - mean(d))[1:(n - j)]) / n, 0)
  expect_equal(r1$other_statistics[["long_run_variance"]],
               g[1] + 2 * sum((1 - seq_len(l) / (l + 1)) * g[-1]), tolerance = 1e-12)
  # trend case
  t <- seq_along(d)
  e <- stats::resid(stats::lm(d ~ t))
  expect_equal(kpss(d, "ct")$other_statistics[["gamma0"]], mean(e^2), tolerance = 1e-10)
})

test_that("large levels do not turn genuine series into constants (major 2)", {
  set.seed(1)
  d <- diff(cumsum(stats::rnorm(200)))
  for (off in c(1e6, 1e12, 1e14, -1e15)) {
    y <- off + d
    # exactly the values stored in y, without the level: the same statistic
    expect_equal(kpss(y)$statistic, kpss(y - off)$statistic, tolerance = 1e-10)
    yt <- y + 0.3 * seq_along(y)
    expect_equal(kpss(yt, "ct")$statistic, kpss(yt - off, "ct")$statistic, tolerance = 1e-10)
  }
  # noise of a few rounding steps on a huge level is still a series; its
  # variance (about 1e570) cannot be represented: NA, statistic unaffected
  off <- 1e300
  y <- off + 8 * off * .Machine$double.eps * stats::rnorm(100)
  r <- kpss(y)
  expect_s3_class(r, "ur_test")
  expect_true(is.na(r$other_statistics[["gamma0"]]))
  expect_true(is.finite(r$statistic))
  # degenerate series are still refused at any level
  for (a in c(1e-200, 1, 1e14, 1e200)) {
    expect_error(kpss(rep(a, 100)), "the series is constant")
    expect_error(kpss(a + 0.5 * a * seq_len(100), "ct"), "exact linear trend")
  }
  expect_error(kpss(rep(0, 100)), "the series is constant")
})
