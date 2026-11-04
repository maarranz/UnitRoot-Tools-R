# Ng-Perron M tests and ERS point-optimal test (Decision D15)

# Independent implementation from Ng and Perron (2001) and Elliott,
# Rothenberg and Stock (1996), written with lm().
gls_indep <- function(y, trend, k) {
  T <- length(y); t <- seq_len(T)
  cbar <- if (trend == "c") -7 else -13.5
  Z <- if (trend == "c") cbind(rep(1, T)) else cbind(1, t)
  qd <- function(x, a) rbind(x[1, , drop = FALSE], x[-1, , drop = FALSE] - a * x[-T, , drop = FALSE])
  a <- 1 + cbar / T
  fa <- stats::lm.fit(qd(Z, a), qd(cbind(y), a)); f1 <- stats::lm.fit(qd(Z, 1), qd(cbind(y), 1))
  yd <- drop(y - Z %*% fa$coefficients)
  dy <- diff(yd); rows <- (k + 1):length(dy)
  X <- cbind(yd[rows]); for (j in seq_len(k)) X <- cbind(X, dy[rows - j])
  f <- stats::lm.fit(X, dy[rows])
  sb <- if (k > 0) sum(f$coefficients[-1]) else 0
  s2 <- (sum(f$residuals^2) / length(rows)) / (1 - sb)^2
  kappa <- sum(yd[-T]^2) / T^2; last <- yd[T]^2 / T
  mza <- (last - s2) / (2 * kappa); msb <- sqrt(kappa / s2)
  mpt <- if (trend == "c") (cbar^2 * kappa - cbar * last) / s2 else (cbar^2 * kappa + (1 - cbar) * last) / s2
  pt <- (sum(fa$residuals^2) - a * sum(f1$residuals^2)) / s2
  c(MZa = mza, MZt = mza * msb, MSB = msb, MPT = mpt, PT = pt)
}

test_that("ng_perron() and ers_pt() follow the published formulas", {
  set.seed(11)
  for (T in c(40, 150)) {
    y <- as.numeric(cumsum(stats::arima.sim(list(ar = 0.4), n = T)))
    for (tr in c("c", "ct")) for (k in c(0L, 3L)) {
      ref <- gls_indep(y, tr, k)
      for (s in c("MZa", "MZt", "MSB", "MPT")) {
        expect_equal(ng_perron(y, tr, lags = k, statistic = s)$statistic, ref[[s]], tolerance = 1e-10)
      }
      expect_equal(ers_pt(y, tr, lags = k)$statistic, ref[["PT"]], tolerance = 1e-10)
    }
  }
})

test_that("GLS tests: lags by Perron-Qu, own critical values and p-values", {
  set.seed(12)
  y <- cumsum(stats::rnorm(200))
  r <- ng_perron(y, "ct")
  sel <- select_lags(y, trend = "ct")
  expect_identical(r$lags, sel$lags)
  expect_identical(r$statistic_name, "MZt")
  expect_named(r$other_statistics, c("MZa", "MSB", "MPT", "s2_AR"))
  expect_identical(r$tail, "lower")
  expect_equal(r$critical_values, rs_critical_values(gls_rs$MZt$ct, 200, r$lags))
  expect_equal(r$p_value, rs_p_value(gls_rs$MZt$ct, r$statistic, 200, r$lags))
  expect_match(r$calibration$source, "own simulation")
  p <- ers_pt(y, "c", lags = 2)
  expect_identical(p$statistic_name, "PT")
  expect_identical(p$test, "ERS-PT")
  for (s in c("MZa", "MZt", "MSB", "MPT", "PT")) for (tr in c("c", "ct")) {
    cv <- rs_critical_values(gls_rs[[s]][[tr]], 100, 2)
    expect_true(cv[[1]] < cv[[2]] && cv[[2]] < cv[[3]])
    pv <- vapply(cv, function(q) rs_p_value(gls_rs[[s]][[tr]], q, 100, 2), numeric(1))
    expect_equal(unname(pv), c(0.01, 0.05, 0.10), tolerance = 0.005)
  }
})

test_that("tabulated calibration: exact on the grid, interpolated in 1/T between", {
  rs <- gls_rs$MZa$c
  j <- which(rs$ks == 2)
  expect_equal(rs_quantiles(rs, 100, 2), unname(rs$table[j, "100", ]))
  q80 <- rs_quantiles(rs, 80, 2); q100 <- rs_quantiles(rs, 100, 2)
  w <- (1 / 90 - 1 / 80) / (1 / 100 - 1 / 80)
  expect_equal(rs_quantiles(rs, 90, 2), (1 - w) * q80 + w * q100)
  expect_equal(rs_quantiles(rs, 1e12, 0), rs$limit, tolerance = 1e-8)
  expect_true(all(diff(rs_quantiles(rs, 37, 5)) > 0))
  # below the design: the smallest simulated T is used, with a warning
  expect_equal(rs_quantiles(rs, 15, 0), rs_quantiles(rs, 20, 0))
  expect_match(rs_design_warnings(rs, 15, 0), "below the smallest sample size")
})

test_that("limits agree with an independent large-T simulation", {
  # 5% quantiles from an independent vectorised simulation, T = 10000, 50000
  # replications (k = 0; spec D15). Ng and Perron's (2001) Table 1 differs in
  # the trend case (MZa -17.3, MZt -2.91), as do ERS's asymptotic DF-GLS
  # values (-2.89); urca's DF-GLS agrees with this package's simulation.
  ref <- list(c = c(MZa = -8.02, MZt = -1.94), ct = c(MZa = -16.69, MZt = -2.853))
  for (tr in c("c", "ct")) for (s in names(ref[[tr]])) {
    cv <- rs_critical_values(gls_rs[[s]][[tr]], 1e9, 0)[["0.05"]]
    expect_equal(cv, ref[[tr]][[s]], tolerance = 0.01, label = paste(s, tr))
  }
  # MPT and PT share their limit
  for (tr in c("c", "ct")) {
    expect_equal(rs_critical_values(gls_rs$MPT[[tr]], 1e9, 0),
                 rs_critical_values(gls_rs$PT[[tr]], 1e9, 0), tolerance = 0.01)
  }
})
