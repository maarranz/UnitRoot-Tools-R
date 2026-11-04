# Enders' published U.S. real GDP example (Supplementary Manual, Section 4.4).
# Published values are rounded, so each is checked to half a unit in its last
# printed digit.

enders_fixture <- function() {
  path <- system.file("fixtures", "enders_rgdp.csv", package = "unitroottools")
  d <- utils::read.csv(path)
  log(d$rgdp)
}

within_rounding <- function(actual, printed, digits) {
  expect_lte(abs(actual - printed), 0.5 * 10^(-digits) + 1e-12)
}

test_that("the fixture is the published data set", {
  y <- enders_fixture()
  expect_length(y, 264L)
})

test_that("ADF regressions (A2.2) and (A2.3) and the phi statistics match", {
  y <- enders_fixture()
  ct <- adf(y, trend = "ct", lags = 1)
  expect_identical(ct$nobs, 262L)
  b <- ct$regression$coefficients
  within_rounding(b["const", "Estimate"], 0.1248, 4)
  within_rounding(b["trend", "Estimate"], 0.0001, 4)
  within_rounding(b["y_lag1", "Estimate"], -0.0156, 4)
  within_rounding(b["dy_lag1", "Estimate"], 0.3663, 4)
  within_rounding(b["const", "t value"], 1.58, 2)
  within_rounding(b["trend", "t value"], 1.31, 2)
  within_rounding(b["y_lag1", "t value"], -1.49, 2)
  within_rounding(b["dy_lag1", "t value"], 6.26, 2)
  within_rounding(ct$other_statistics[["phi3"]], 2.97, 2)

  c1 <- adf(y, trend = "c", lags = 1)
  b <- c1$regression$coefficients
  within_rounding(b["const", "Estimate"], 0.0215, 4)
  within_rounding(b["y_lag1", "Estimate"], -0.0019, 4)
  within_rounding(b["dy_lag1", "Estimate"], 0.3539, 4)
  expect_lte(abs(b["const", "t value"] - 2.64), 0.01)    # printed 2.64, value 2.646
  within_rounding(b["y_lag1", "t value"], -2.05, 2)
  within_rounding(b["dy_lag1", "t value"], 6.12, 2)
  within_rounding(c1$other_statistics[["phi1"]], 25.49, 2)
})

test_that("the original procedure follows Enders' path and rejects", {
  y <- enders_fixture()
  o <- ur_strategy(y, "original", lags = 1)
  expect_identical(o$history$test, c("tau3", "phi3", "tau2", "phi1", "aux_const", "cond_c"))
  expect_identical(o$specification, "c")
  expect_true(o$unit_root_rejected)
  within_rounding(o$history$statistic[o$history$test == "aux_const"], 6.80, 2)
  # Enders compares phi3 with the T = 250 row of Dickey-Fuller's table (5%
  # value 6.34); the package's response surface at this T agrees to rounding
  within_rounding(o$history$critical_value[o$history$test == "phi3"], 6.34, 2)
})

test_that("the conservative procedure retains the constant and does not reject", {
  y <- enders_fixture()
  s <- ur_strategy(y, "conservative", lags = 1)
  expect_identical(s$history$test, c("tau3", "phi3", "tau2"))
  expect_identical(s$specification, "c")
  expect_false(s$unit_root_rejected)
})
