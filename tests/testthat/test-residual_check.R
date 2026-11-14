# Residual diagnostics (spec: Decision D19)

ar1_rw <- function(n = 240, phi = 0.6, freq = 12) {
  set.seed(42)
  stats::ts(cumsum(stats::arima.sim(list(ar = phi), n = n)), start = 2000, frequency = freq)
}

test_that("Breusch-Godfrey agrees with lmtest::bgtest (F and LM forms)", {
  skip_if_not_installed("lmtest")
  y <- ar1_rw()
  for (spec in list(list(trend = "c", k = 0L), list(trend = "ct", k = 2L), list(trend = "n", k = 1L))) {
    r <- adf(y, trend = spec$trend, lags = spec$k)
    X <- r$regression$X
    resp <- r$regression$residuals + drop(X %*% r$regression$coefficients[, "Estimate"])
    fit <- stats::lm(resp ~ X - 1)
    rc <- residual_check(r)
    for (m in rc$breusch_godfrey$order) {
      ref_f <- lmtest::bgtest(fit, order = m, type = "F")
      ref_c <- lmtest::bgtest(fit, order = m, type = "Chisq")
      row <- rc$breusch_godfrey[rc$breusch_godfrey$order == m, ]
      expect_equal(row$f, unname(ref_f$statistic), tolerance = 1e-8)
      expect_equal(c(row$df1, row$df2), unname(ref_f$parameter), tolerance = 1e-12)
      expect_equal(row$p_value_f, ref_f$p.value, tolerance = 1e-8)
      expect_equal(row$lm, unname(ref_c$statistic), tolerance = 1e-8)
      expect_equal(row$p_value_lm, ref_c$p.value, tolerance = 1e-8)
    }
  }
})

test_that("correlogram agrees with stats::acf, pacf and Box.test", {
  y <- ar1_rw()
  r <- adf(y, trend = "c", lags = 2)
  rc <- residual_check(r)
  e <- r$regression$residuals
  expect_equal(rc$correlogram$ac, drop(stats::acf(e, lag.max = rc$lags, plot = FALSE)$acf)[-1])
  expect_equal(rc$correlogram$pac, drop(stats::pacf(e, lag.max = rc$lags, plot = FALSE)$acf))
  for (h in c(3, 12, 39)) {
    b <- stats::Box.test(e, lag = h, type = "Ljung-Box", fitdf = 2)
    expect_equal(rc$correlogram$q[h], unname(b$statistic), tolerance = 1e-10)
    expect_equal(rc$correlogram$p_value[h], b$p.value, tolerance = 1e-10)
  }
  expect_true(all(is.na(rc$correlogram$p_value[1:2])))   # h <= k: no df
})

test_that("lag grid follows the frequency: 3s + 3, capped at N/4", {
  y <- ar1_rw()
  rc <- residual_check(adf(y, lags = 1))
  expect_identical(rc$frequency, 12L)
  expect_identical(rc$lags, 39L)
  expect_identical(rc$breusch_godfrey$order, c(12L, 39L))
  q <- ar1_rw(freq = 4)
  rq <- residual_check(adf(q, lags = 1))
  expect_identical(c(rq$lags, rq$breusch_godfrey$order), c(15L, 4L, 15L))
  # plain vector: annual with a note; frequency argument overrides
  v <- as.numeric(y)
  rv <- residual_check(adf(v, lags = 1))
  expect_identical(c(rv$frequency, rv$lags), c(1L, 6L))
  expect_match(rv$notes, "frequency unknown", all = FALSE)
  expect_identical(residual_check(adf(v, lags = 1), frequency = 12)$lags, 39L)
  # short sample: capped at N/4
  short <- stats::ts(as.numeric(y)[1:80], frequency = 12)
  rs <- residual_check(adf(short, lags = 1))
  expect_identical(rs$lags, 78L %/% 4L)
  expect_match(rs$notes, "capped at N/4", all = FALSE)
  expect_true(all(rs$breusch_godfrey$order <= rs$lags))
  # explicit choices
  rx <- residual_check(adf(y, lags = 1), lags = 24, bg_orders = c(1, 6))
  expect_identical(c(rx$lags, rx$breusch_godfrey$order), c(24L, 1L, 6L))
})

test_that("missing lags are detected and enough lags pass", {
  y <- ar1_rw(phi = 0.7)
  bad <- residual_check(adf(y, "c", lags = 0))
  good <- residual_check(adf(y, "c", lags = 1))
  expect_false(bad$white_noise)
  expect_true(good$white_noise)
  expect_match(capture.output(print(bad)), "try a larger `lags`", all = FALSE)
  expect_match(capture.output(print(good)), "No evidence of residual autocorrelation", all = FALSE)
})

test_that("GLS tests are checked on their own regression; PP and KPSS are refused", {
  y <- ar1_rw()
  for (r in list(dfgls(y, "c", lags = 1), ng_perron(y, "ct", lags = 1), ers_pt(y, "c", lags = 1))) {
    rc <- residual_check(r)
    expect_s3_class(rc, "ur_residual_check")
    expect_identical(rc$nobs, r$nobs)
  }
  expect_error(residual_check(pp(y)), "by design")
  expect_error(residual_check(kpss(y)), "by design")
  old <- adf(y, lags = 1); old$regression$residuals <- NULL
  expect_error(residual_check(old), "does not contain the residuals")
})

test_that("plot draws without error", {
  rc <- residual_check(adf(ar1_rw(), lags = 1))
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f); on.exit(unlink(f))
  expect_invisible(plot(rc))
  grDevices::dev.off()
})
