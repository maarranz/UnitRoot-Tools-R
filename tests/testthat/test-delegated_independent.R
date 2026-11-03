# Phase 2B audit: the statistics computed by urca are checked against
# independent implementations written from the published formulas.

bartlett_lrv <- function(u, l) {
  n <- length(u); g0 <- sum(u^2) / n; s <- g0
  for (j in seq_len(l)) s <- s + 2 * (1 - j / (l + 1)) * sum(u[(j + 1):n] * u[1:(n - j)]) / n
  c(g0 = g0, lrv = s)
}

# KPSS (1992): eta = T^{-2} sum S_t^2 / lambda^2
kpss_indep <- function(y, trend, l) {
  T <- length(y); t <- seq_len(T)
  e <- if (trend == "c") y - mean(y) else stats::lm.fit(cbind(1, t), y)$residuals
  sum(cumsum(e)^2) / (T^2 * bartlett_lrv(e, l)[["lrv"]])
}

# Elliott, Rothenberg and Stock (1996): GLS detrending with cbar = -7 / -13.5,
# then the t ratio of yd_{t-1} in a regression without deterministic terms.
dfgls_indep <- function(y, trend, k) {
  T <- length(y); t <- seq_len(T)
  a <- 1 + (if (trend == "c") -7 else -13.5) / T
  Z <- if (trend == "c") cbind(rep(1, T)) else cbind(1, t)
  qd <- function(x) rbind(x[1, , drop = FALSE], x[-1, , drop = FALSE] - a * x[-T, , drop = FALSE])
  yd <- drop(y - Z %*% qr.coef(qr(qd(Z)), qd(cbind(y))))
  dy <- diff(yd); rows <- (k + 1):length(dy)
  X <- cbind(yd[rows])
  for (j in seq_len(k)) X <- cbind(X, dy[rows - j])
  f <- stats::lm.fit(X, dy[rows])
  s2 <- sum(f$residuals^2) / (length(rows) - ncol(X))
  unname(f$coefficients[1] / sqrt(s2 * chol2inv(qr.R(qr(X)))[1, 1]))
}

# Phillips-Perron (1988) Z_tau with a constant, written two ways: with the
# published sums over y_{t-1}, and with urca's sums over y_t.
pp_c <- function(x, l, lagged) {
  n <- length(x) - 1; y <- x[-1]; yl <- x[-(n + 1)]
  f <- summary(stats::lm(y ~ yl))
  tstat <- (stats::coef(f)[2, 1] - 1) / stats::coef(f)[2, 2]
  v <- bartlett_lrv(stats::residuals(f), l)
  z <- if (lagged) yl else y
  sqrt(v[["g0"]] / v[["lrv"]]) * tstat -
    0.5 * (v[["lrv"]] - v[["g0"]]) / sqrt(v[["lrv"]]) / sqrt(sum((z - mean(z))^2) / n^2)
}

test_that("kpss() and dfgls() statistics match the published formulas", {
  set.seed(42)
  for (T in c(30, 120, 264)) {
    y <- cumsum(0.1 + stats::arima.sim(list(ar = 0.5), n = T))
    y <- as.numeric(y)
    for (tr in c("c", "ct")) {
      k <- kpss(y, trend = tr)
      expect_equal(k$statistic, kpss_indep(y, tr, k$lags), tolerance = 1e-10)
      for (p in c(0, 1, 4)) {
        expect_equal(dfgls(y, trend = tr, lags = p)$statistic, dfgls_indep(y, tr, p),
                     tolerance = 1e-10)
      }
    }
  }
})

# Phillips-Perron Z_tau written independently (Hamilton 1994, eq. 17.6.8):
# levels regression y_t on (1, [t - n/2], y_{t-1}).
pp_indep <- function(y, trend, l) {
  T <- length(y); yt <- y[-1]; yl <- y[-T]; n <- T - 1
  X <- if (trend == "c") cbind(1, yl) else cbind(1, seq_len(n) - n / 2, yl)
  f <- stats::lm.fit(X, yt); u <- f$residuals; k <- ncol(X)
  s2 <- sum(u^2) / (n - k)
  se_rho <- sqrt(s2 * chol2inv(qr.R(qr(X)))[k, k])
  t_rho <- (f$coefficients[k] - 1) / se_rho
  v <- bartlett_lrv(u, l)
  unname(sqrt(v[["g0"]] / v[["lrv"]]) * t_rho -
           0.5 * (v[["lrv"]] - v[["g0"]]) / sqrt(v[["lrv"]]) * n * se_rho / sqrt(s2))
}

test_that("pp() follows the published Phillips-Perron formula (D12)", {
  set.seed(43)
  for (T in c(30, 120, 264)) {
    y <- as.numeric(cumsum(0.1 + stats::arima.sim(list(ma = -0.4), n = T)))
    for (tr in c("c", "ct")) for (bw in list("short", "long", 0L, 5L)) {
      p <- pp(y, trend = tr, bandwidth = bw)
      expect_equal(p$statistic, pp_indep(y, tr, p$lags), tolerance = 1e-10)
    }
  }
  # urca's version (sums over y_t) differs visibly at T = 30
  y <- cumsum(stats::rnorm(30))
  p <- pp(y, trend = "c")
  expect_equal(p$statistic, pp_c(y, p$lags, lagged = TRUE), tolerance = 1e-10)
  expect_gt(abs(p$statistic - pp_c(y, p$lags, lagged = FALSE)), 1e-3)
})
