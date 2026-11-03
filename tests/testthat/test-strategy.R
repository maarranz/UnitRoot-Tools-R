# ---- decision logic: every branch, driven by constructed decisions ----------

rej <- function(...) {
  r <- c(tau3 = FALSE, phi3 = FALSE, aux_trend = FALSE, cond_ct = FALSE,
         tau2 = FALSE, phi1 = FALSE, aux_const = FALSE, cond_c = FALSE,
         tau1 = FALSE)
  set <- c(...)
  r[names(set)] <- set
  r
}
visited <- function(p) p$steps$test

test_that("original: tau3 rejects -> stop with trend", {
  p <- strategy_path(rej(tau3 = TRUE), "original")
  expect_identical(visited(p), "tau3")
  expect_identical(p$specification, "ct")
  expect_true(p$unit_root_rejected)
})

test_that("original: phi3 does not reject -> constant-only stage", {
  p <- strategy_path(rej(tau2 = TRUE), "original")
  expect_identical(visited(p), c("tau3", "phi3", "tau2"))
  expect_identical(p$terminal, "tau2")
  expect_identical(p$specification, "c")
  expect_true(p$unit_root_rejected)
})

test_that("original: phi3 rejects, trend not significant -> constant-only stage", {
  p <- strategy_path(rej(phi3 = TRUE), "original")
  expect_identical(visited(p)[1:4], c("tau3", "phi3", "aux_trend", "tau2"))
})

test_that("original: trend significant -> conditional t test of b (both outcomes)", {
  p <- strategy_path(rej(phi3 = TRUE, aux_trend = TRUE, cond_ct = TRUE), "original")
  expect_identical(visited(p), c("tau3", "phi3", "aux_trend", "cond_ct"))
  expect_identical(p$terminal, "cond_ct")
  expect_identical(p$specification, "ct")
  expect_true(p$unit_root_rejected)
  p2 <- strategy_path(rej(phi3 = TRUE, aux_trend = TRUE), "original")
  expect_false(p2$unit_root_rejected)
  expect_match(p2$conclusion, "quadratic trend")
})

test_that("original: constant stage branches", {
  # phi1 does not reject -> tau1
  p <- strategy_path(rej(), "original")
  expect_identical(visited(p), c("tau3", "phi3", "tau2", "phi1", "tau1"))
  expect_identical(p$specification, "n")
  expect_false(p$unit_root_rejected)
  expect_match(p$conclusion, "without drift")
  # phi1 rejects, drift not significant -> tau1
  p <- strategy_path(rej(phi1 = TRUE, tau1 = TRUE), "original")
  expect_identical(visited(p), c("tau3", "phi3", "tau2", "phi1", "aux_const", "tau1"))
  expect_true(p$unit_root_rejected)
  # phi1 rejects, drift significant -> conditional t (Enders' GDP path)
  p <- strategy_path(rej(phi1 = TRUE, aux_const = TRUE, cond_c = TRUE), "original")
  expect_identical(visited(p), c("tau3", "phi3", "tau2", "phi1", "aux_const", "cond_c"))
  expect_identical(p$specification, "c")
  expect_true(p$unit_root_rejected)
  p <- strategy_path(rej(phi1 = TRUE, aux_const = TRUE), "original")
  expect_identical(p$conclusion, "unit root with drift")
})

test_that("conservative: constant retained when tau2 does not reject", {
  # Same decisions as Enders' GDP path: the original rejects, the conservative does not
  d <- rej(phi1 = TRUE, aux_const = TRUE, cond_c = TRUE)
  o <- strategy_path(d, "original")
  c1 <- strategy_path(d, "conservative")
  expect_true(o$unit_root_rejected)
  expect_identical(visited(c1), c("tau3", "phi3", "tau2"))
  expect_identical(c1$terminal, "tau2")
  expect_identical(c1$specification, "c")
  expect_false(c1$unit_root_rejected)
  # Stage T is identical to the original
  expect_identical(strategy_path(rej(tau3 = TRUE), "conservative")$terminal, "tau3")
  expect_identical(visited(strategy_path(rej(phi3 = TRUE, aux_trend = TRUE), "conservative")),
                   c("tau3", "phi3", "aux_trend", "cond_ct"))
})

test_that("conservative with allow_no_constant follows the original constant stage", {
  d <- rej(phi1 = TRUE, aux_const = TRUE, cond_c = TRUE)
  expect_identical(visited(strategy_path(d, "conservative", allow_no_constant = TRUE)),
                   visited(strategy_path(d, "original")))
  d2 <- rej()
  expect_identical(strategy_path(d2, "conservative", allow_no_constant = TRUE)$terminal, "tau1")
})

test_that("hamilton: one terminal test chosen by the declaration", {
  p <- strategy_path(rej(tau3 = TRUE), "hamilton", trending = TRUE)
  expect_identical(visited(p), "tau3")
  expect_identical(p$specification, "ct")
  expect_true(p$unit_root_rejected)
  p <- strategy_path(rej(tau3 = TRUE), "hamilton", trending = FALSE)
  expect_identical(visited(p), "tau2")
  expect_false(p$unit_root_rejected)
})

# ---- statistics on real computations -----------------------------------------

sim_trend <- function(n = 200, seed = 4) {
  set.seed(seed)
  cumsum(0.3 + stats::arima.sim(list(ar = 0.3), n = n))
}

test_that("strategy statistics equal adf() on the same sample and lag order", {
  y <- sim_trend()
  st <- strategy_statistics(y, k = 2, first = 4, level = 0.05)
  a_ct <- adf(y, trend = "ct", lags = 2)
  a_c <- adf(y, trend = "c", lags = 2)
  a_n <- adf(y, trend = "n", lags = 2)
  expect_equal(st$tests$tau3$statistic, a_ct$statistic, tolerance = 1e-12)
  expect_equal(st$tests$phi3$statistic, a_ct$other_statistics[["phi3"]], tolerance = 1e-12)
  expect_equal(st$tests$tau2$statistic, a_c$statistic, tolerance = 1e-12)
  expect_equal(st$tests$phi1$statistic, a_c$other_statistics[["phi1"]], tolerance = 1e-12)
  expect_equal(st$tests$tau1$statistic, a_n$statistic, tolerance = 1e-12)
  expect_equal(st$tests$tau3$critical_value, a_ct$critical_values[["0.05"]])
})

test_that("auxiliary and conditional tests match lm() with Student t", {
  y <- sim_trend(seed = 9)
  k <- 1; first <- 3
  st <- strategy_statistics(y, k, first, level = 0.05)
  t <- first:length(y); dy <- diff(y)
  dat <- data.frame(dy = dy[t - 1], tt = t, ylag = y[t - 1], l1 = dy[t - 2])
  aux_ct <- summary(stats::lm(dy ~ tt + l1, dat))$coefficients
  expect_equal(st$tests$aux_trend$statistic, aux_ct["tt", 3], tolerance = 1e-10)
  aux_c <- summary(stats::lm(dy ~ l1, dat))$coefficients
  expect_equal(st$tests$aux_const$statistic, aux_c["(Intercept)", 3], tolerance = 1e-10)
  expect_equal(st$tests$aux_trend$critical_value, stats::qt(0.975, nrow(dat) - 3))
  full_c <- stats::lm(dy ~ ylag + l1, dat)
  expect_equal(st$tests$cond_c$statistic,
               summary(full_c)$coefficients["ylag", 3], tolerance = 1e-10)
  expect_equal(st$tests$cond_c$critical_value, stats::qt(0.05, full_c$df.residual))
})

test_that("ur_strategy: one lag order and one sample for every step", {
  y <- sim_trend()
  s <- ur_strategy(y, "original")
  expect_s3_class(s, "ur_strategy")
  expect_identical(s$lag_selection$selected_under, "ct")
  expect_identical(s$lags, select_lags(y, trend = "ct")$lags)
  expect_identical(s$nobs, length(y) - s$lags - 1L)
  expect_identical(s$lag_selection$final_sample, "maximal")
  sc <- ur_strategy(y, "original", final_sample = "common")
  expect_identical(sc$nobs, select_lags(y, trend = "ct")$nobs)
  expect_identical(sc$terminal$nobs, sc$nobs)
  expect_s3_class(s$terminal, "ur_test")
  expect_identical(s$terminal$nobs, s$nobs)
  expect_identical(s$history$step, seq_len(nrow(s$history)))
  expect_true(all(s$history$rejection_type[!s$history$reject] == "none"))
  expect_identical(as.data.frame(s), s$history)
})

test_that("ur_strategy: rejection types separate unit-root and path decisions", {
  y <- sim_trend()
  s <- ur_strategy(y, "original", lags = 1)
  h <- s$history
  types <- unname(rejection_types[h$test])
  expect_identical(h$rejection_type, ifelse(h$reject, types, "none"))
})

test_that("ur_strategy: the terminal test agrees with the history", {
  for (seed in 1:6) {
    y <- sim_trend(seed = seed)
    for (strategy in c("original", "conservative")) {
      s <- ur_strategy(y, strategy, lags = 1)
      last <- s$history[nrow(s$history), ]
      expect_identical(s$unit_root_rejected, last$reject)
      expect_equal(s$terminal$statistic, last$statistic)
      expect_identical(s$terminal$reject, last$reject)
    }
  }
})

test_that("ur_strategy: conditional terminal test uses Student t with a p-value", {
  y <- sim_trend()
  st <- strategy_statistics(y, 1, 3, 0.05)
  tt <- terminal_test("cond_c", st, NULL, 0.05, NULL)
  expect_identical(tt$calibration$method, "conventional_t")
  expect_equal(tt$p_value, stats::pt(tt$statistic, st$tests$cond_c$df))
})

test_that("hamilton requires a declaration; misuse is reported", {
  y <- sim_trend()
  expect_error(ur_strategy(y, "hamilton"), "There is no default")
  expect_error(ur_strategy(y, "hamilton", trending = NA), "There is no default")
  h <- ur_strategy(y, "hamilton", trending = TRUE)
  expect_identical(h$history$test, "tau3")
  expect_identical(h$lag_selection$selected_under, "ct")
  expect_identical(ur_strategy(y, "hamilton", trending = FALSE)$lag_selection$selected_under, "c")
  expect_warning(ur_strategy(y, "original", trending = TRUE), "only used by")
  expect_warning(ur_strategy(y, "original", allow_no_constant = TRUE), "only used by")
})

test_that("argument checks", {
  y <- sim_trend()
  expect_error(ur_strategy(y, level = 0.025), "one of 0.01, 0.05, 0.1")
  expect_error(ur_strategy(y[1:30], lags = 8), "at least 25 effective")
  expect_error(ur_strategy(y, lags = 1, max_lags = 4), "either `lags`")
})

test_that("print and summary describe the path", {
  y <- sim_trend()
  s <- ur_strategy(y, "original")
  out <- capture.output(print(s))
  expect_match(out, "Specification strategy: Original", all = FALSE)
  expect_match(out, "Conclusion:", all = FALSE)
  expect_match(out, "no p-value or size is reported", all = FALSE)
  sm <- capture.output(print(summary(s)))
  expect_match(sm, "H0: b = 0: unit root \\(ct regression\\)", all = FALSE)
  expect_match(sm, "Terminal test", all = FALSE)
})
