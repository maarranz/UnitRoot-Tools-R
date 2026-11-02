sim_ar_rw <- function(n = 200, seed = 11) {
  set.seed(seed)
  cumsum(0.1 + stats::arima.sim(list(ar = c(0.5, -0.2)), n = n))
}

test_that("all orders are compared on the common sample", {
  y <- sim_ar_rw()
  s <- select_lags(y, trend = "ct", max_lags = 8)
  expect_s3_class(s, "ur_lags")
  expect_identical(s$kmax, 8L)
  expect_identical(s$nobs, 200L - 8L - 1L)
  expect_identical(as.integer(s$sample), c(10L, 200L))
  expect_identical(s$table$k, 0:8)
  expect_named(s$table, c("k", "AIC", "BIC", "MAIC", "MBIC"))
})

test_that("default kmax is Schwert's rule", {
  y <- sim_ar_rw(n = 200)
  expect_identical(select_lags(y)$kmax, as.integer(floor(12 * (200 / 100)^0.25)))
})

test_that("AIC and BIC match lm() fits on the common sample", {
  y <- sim_ar_rw()
  kmax <- 6
  s <- select_lags(y, trend = "c", max_lags = kmax, criterion = "AIC")
  n <- length(y); t <- (kmax + 2):n; N <- length(t)
  dy <- diff(y)
  for (k in 0:kmax) {
    lagm <- if (k > 0) sapply(seq_len(k), function(j) dy[t - 1 - j]) else NULL
    X <- cbind(ylag = y[t - 1], lagm)
    rss <- sum(stats::resid(stats::lm(dy[t - 1] ~ X))^2)
    expect_equal(s$table$AIC[k + 1], log(rss / N) + 2 * k / N, tolerance = 1e-10)
    expect_equal(s$table$BIC[k + 1], log(rss / N) + k * log(N) / N, tolerance = 1e-10)
  }
  expect_identical(s$lags, s$table$k[which.min(s$table$AIC)])
})

test_that("MAIC and MBIC match an independent Ng-Perron computation", {
  y <- sim_ar_rw(seed = 5)
  kmax <- 6
  n <- length(y); tt <- seq_len(n)
  yt <- stats::resid(stats::lm(y ~ tt))               # OLS detrending ("ct")
  t <- (kmax + 2):n; N <- length(t)
  dyt <- diff(yt)
  s <- select_lags(y, trend = "ct", max_lags = kmax)
  for (k in 0:kmax) {
    lagm <- if (k > 0) sapply(seq_len(k), function(j) dyt[t - 1 - j]) else NULL
    X <- cbind(ylag = yt[t - 1], lagm)
    f <- stats::lm(dyt[t - 1] ~ X - 1)
    s2 <- sum(stats::resid(f)^2) / N
    tau <- stats::coef(f)[1]^2 * sum(yt[t - 1]^2) / s2
    expect_equal(s$table$MAIC[k + 1], unname(log(s2) + 2 * (tau + k) / N), tolerance = 1e-10)
    expect_equal(s$table$MBIC[k + 1], unname(log(s2) + log(N) * (tau + k) / N), tolerance = 1e-10)
  }
})

test_that("criteria, input checks and printing", {
  y <- sim_ar_rw()
  expect_identical(select_lags(y, criterion = "BIC")$criterion, "BIC")
  expect_error(select_lags(y, criterion = "HQ"), "should be one of")
  expect_error(select_lags(y[1:10], max_lags = 8), "too short")
  expect_error(select_lags(y, max_lags = -1), "non-negative")
  expect_identical(select_lags(y, trend = "n")$detrending, "none")
  out <- capture.output(print(select_lags(y, max_lags = 4)))
  expect_match(out, "Common sample:  N = 195", all = FALSE)
  expect_match(out, "Selected by MAIC: k = ", all = FALSE)
})
