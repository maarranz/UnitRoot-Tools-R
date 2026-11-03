# Long-run variance module (Decision D13)

test_that("kernel weights", {
  x <- c(0, 0.25, 0.5, 0.75, 1, 1.5)
  expect_equal(kernel_weight(x, "bartlett"), c(1, 0.75, 0.5, 0.25, 0, 0))
  expect_equal(kernel_weight(x, "parzen"),
               c(1, 1 - 6 * 0.25^2 + 6 * 0.25^3, 2 * 0.5^3, 2 * 0.25^3, 0, 0))
  z <- 6 * pi * 0.5 / 5
  expect_equal(kernel_weight(0.5, "qs"), 25 / (12 * pi^2 * 0.25) * (sin(z) / z - cos(z)))
  expect_equal(kernel_weight(0, "qs"), 1)
})

test_that("Schwert bandwidths reproduce the Bartlett weights 1 - j/(l+1)", {
  set.seed(1)
  u <- as.numeric(stats::arima.sim(list(ar = 0.5), n = 200)); u <- u - mean(u)
  bw <- resolve_bandwidth("short", u, "bartlett")
  expect_identical(bw$lags, 4L)
  g <- vapply(0:4, function(j) sum(u[(j + 1):200] * u[1:(200 - j)]) / 200, 0)
  expect_equal(long_run_variance(u, "bartlett", bw$b)$lrv,
               g[1] + 2 * sum((1 - (1:4) / 5) * g[-1]))
  expect_equal(long_run_variance(u, "bartlett", 1)$lrv, g[1])     # 0 autocovariances
})

test_that("long-run variances and Newey-West bandwidths agree with sandwich", {
  skip_if_not_installed("sandwich")
  set.seed(2)
  for (r in 1:8) {
    n <- c(60, 150, 400)[(r %% 3) + 1]
    u <- as.numeric(stats::arima.sim(list(ar = c(-0.4, 0.3, 0.8)[(r %% 3) + 1]), n = n))
    u <- u - mean(u)
    m <- stats::lm(u ~ 1)
    for (k in lrv_kernels) {
      sk <- c(bartlett = "Bartlett", parzen = "Parzen", qs = "Quadratic Spectral")[[k]]
      expect_equal(newey_west_bandwidth(u, k),
                   sandwich::bwNeweyWest(m, kernel = sk, prewhite = 0, adjust = FALSE),
                   tolerance = 1e-10)
      b <- 1 + 7 * r / 8
      expect_equal(long_run_variance(u, k, b)$lrv,
                   drop(sandwich::kernHAC(m, kernel = sk, bw = b, prewhite = 0,
                                          adjust = FALSE, sandwich = FALSE)),
                   tolerance = 1e-10)
      # Andrews' AR(1) is fitted without an intercept (sandwich demeans first)
      expect_equal(andrews_bandwidth(u, k),
                   sandwich::bwAndrews(m, kernel = sk, approx = "AR(1)", prewhite = 0,
                                       adjust = FALSE), tolerance = 0.02)
    }
  }
})

test_that("KPSS limiting distributions: exact eigenvalues and quantiles", {
  expect_equal(kpss_limit$mu[1:5], 1 / ((1:5) * pi)^2, tolerance = 1e-6)
  expect_equal(sum(kpss_limit$mu), 1 / 6, tolerance = 1e-6)
  expect_equal(sum(kpss_limit$tau), 1 / 15, tolerance = 1e-6)
  for (tr in c("c", "ct")) {
    cv <- kpss_critical_values(tr)
    lam <- kpss_limit[[c(c = "mu", ct = "tau")[[tr]]]]
    p <- vapply(cv, kpss_upper_prob, numeric(1), lam = lam)
    expect_equal(unname(p), c(0.10, 0.05, 0.025, 0.01), tolerance = 1e-6)
  }
  expect_equal(unname(kpss_critical_values("c")), c(0.3473, 0.4614, 0.5806, 0.7435), tolerance = 2e-4)
  expect_equal(unname(kpss_critical_values("ct")), c(0.1192, 0.1479, 0.1775, 0.2177), tolerance = 2e-4)
})

test_that("kpss() and pp() accept every kernel and automatic bandwidth", {
  set.seed(3)
  y <- cumsum(stats::rnorm(120))
  for (k in lrv_kernels) for (b in list("short", "long", "andrews", "newey-west", 3L)) {
    r <- kpss(y, trend = "ct", kernel = k, bandwidth = b)
    expect_s3_class(r, "ur_test")
    expect_match(r$bandwidth_rule, c(bartlett = "Bartlett", parzen = "Parzen", qs = "Quadratic Spectral")[[k]])
    p <- pp(y, trend = "c", kernel = k, bandwidth = b)
    expect_s3_class(p, "ur_test")
  }
  # Bartlett with a fixed number of autocovariances: unchanged definitions
  expect_equal(kpss(y, bandwidth = 4)$statistic, kpss(y)$statistic)
})

test_that("pp(): Z_alpha follows Hamilton (1994, eq. 17.6.7)", {
  set.seed(4)
  y <- cumsum(stats::rnorm(90))
  p <- pp(y, trend = "c", bandwidth = 3)
  n <- 89; yt <- y[-1]; yl <- y[-90]
  f <- stats::lm(yt ~ yl); u <- stats::residuals(f)
  s2 <- sum(u^2) / (n - 2); se <- summary(f)$coefficients["yl", 2]
  g <- vapply(0:3, function(j) sum(u[(j + 1):n] * u[1:(n - j)]) / n, 0)
  lam2 <- g[1] + 2 * sum((1 - (1:3) / 4) * g[-1])
  za <- n * (stats::coef(f)[["yl"]] - 1) - 0.5 * (n^2 * se^2 / s2) * (lam2 - g[1])
  expect_equal(p$other_statistics[["Z_alpha"]], za, tolerance = 1e-10)
})
