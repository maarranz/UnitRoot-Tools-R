# DF-GLS critical values and p-values from the package's own response
# surfaces (Decision D11).

test_that("critical values: ordered, consistent with p-values, right limits", {
  for (tr in c("c", "ct")) for (T in c(25, 40, 100, 264, 1000)) for (k in c(0, 1, 5, 12, 21)) {
    if (k > floor(12 * (T / 100)^(1 / 4))) next          # inside the design only
    cv <- dfgls_critical_values(tr, T, k)
    expect_named(cv, c("0.01", "0.05", "0.1"))
    expect_true(cv[[1]] < cv[[2]] && cv[[2]] < cv[[3]])
    p <- vapply(cv, dfgls_p_value, numeric(1), trend = tr, T = T, k = k)
    expect_equal(unname(p), c(0.01, 0.05, 0.10), tolerance = 0.005)
  }
  # asymptotics: with a constant, the Dickey-Fuller distribution without
  # deterministic terms (5%: -1.941); both quantiles flat in k as T grows
  inf_c <- dfgls_rs$c$intercept[which.min(abs(dfgls_rs$c$probs - 0.05))]
  expect_equal(inf_c, -1.941, tolerance = 0.01)
  expect_lt(abs(dfgls_critical_values("ct", 1e6, 0)[[2]] -
                  dfgls_critical_values("ct", 1e6, 12)[[2]]), 1e-3)
})

test_that("finite-sample critical values are more conservative than the asymptotic ones", {
  expect_lt(dfgls_critical_values("c", 50, 0)[["0.05"]], -2.2)
  expect_lt(dfgls_critical_values("c", 100, 0)[["0.05"]], -2.1)
  expect_lt(dfgls_critical_values("ct", 50, 0)[["0.05"]], -3.15)
})

test_that("p-values are monotone in the statistic and stay in [0, 1]", {
  s <- seq(-6, 3, by = 0.25)
  p <- vapply(s, dfgls_p_value, numeric(1), trend = "c", T = 120, k = 2)
  expect_true(all(diff(p) >= -1e-8))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("dfgls() reports own critical values, p-value and design warnings", {
  set.seed(3)
  y <- cumsum(stats::rnorm(150))
  r <- dfgls(y, trend = "ct", lags = 2)
  expect_equal(r$critical_values, dfgls_critical_values("ct", 150, 2))
  expect_equal(r$p_value, dfgls_p_value(r$statistic, "ct", 150, 2))
  expect_identical(r$reject, r$statistic < r$critical_values[["0.05"]])
  expect_match(r$calibration$source, "own simulation")
  expect_length(r$warnings, 0L)
  expect_match(dfgls(y[1:22], trend = "c", lags = 1)$warnings, "below the smallest", all = FALSE)
  # outside the design the surface is not extrapolated
  expect_equal(dfgls_critical_values("c", 40, 12), dfgls_critical_values("c", 100, 12))
  expect_match(dfgls(y, trend = "c", lags = 25)$warnings, "exceeds the largest", all = FALSE)
})
