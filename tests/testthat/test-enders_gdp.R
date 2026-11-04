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

# The shared fixture file (audit N4): settings and expected values are read
# from inst/fixtures/enders_rgdp.json, so R and Python check the same numbers.
test_that("results match the shared Enders fixture (enders_rgdp.json)", {
  skip_if_not_installed("jsonlite")
  fx <- jsonlite::read_json(system.file("fixtures", "enders_rgdp.json", package = "unitroottools"))
  ex <- fx$expected
  # data identity: the CSV whose SHA-256 the fixture records
  # (9a857973...d74d9; base R has only MD5, recorded here for the same file)
  expect_identical(fx$source$csv_sha256,
                   "9a857973fcf968ac62bcce3e9ace55d52af993c6003169178706e39daadd74d9")
  expect_identical(unname(tools::md5sum(system.file("fixtures", fx$source$data_file,
                                                    package = "unitroottools"))),
                   "baff50bff996f44e2ccd0401828ec848")
  y <- enders_fixture()
  k <- fx$settings$lags
  ct <- adf(y, trend = "ct", lags = k)
  c1 <- adf(y, trend = "c", lags = k)
  expect_identical(ct$nobs, as.integer(fx$settings$nobs))
  check_block <- function(fit, coefs, column, digits) {
    for (nm in names(coefs)) {
      printed <- coefs[[nm]]
      actual <- fit$regression$coefficients[nm, column]
      # the documented exception: the (A2.3) intercept t value is truncated
      if (identical(column, "t value") && nm == "const" && identical(fit$trend, "c")) {
        expect_lte(abs(actual - printed), 0.01)
      } else {
        within_rounding(actual, printed, digits)
      }
    }
  }
  check_block(ct, ex$ct_coefficients, "Estimate", 4)
  check_block(ct, ex$ct_t_values, "t value", 2)
  check_block(c1, ex$c_coefficients, "Estimate", 4)
  check_block(c1, ex$c_t_values, "t value", 2)
  within_rounding(ct$other_statistics[["phi3"]], ex$phi3, 2)
  within_rounding(c1$other_statistics[["phi1"]], ex$phi1, 2)
  for (strategy in c("original", "conservative")) {
    s <- ur_strategy(y, strategy, lags = k, level = fx$settings$level)
    expect_identical(s$history$test, unlist(ex[[strategy]]$tests))
    expect_identical(s$unit_root_rejected, ex[[strategy]]$unit_root_rejected)
    expect_identical(s$specification, ex[[strategy]]$specification)
  }
  o <- ur_strategy(y, "original", lags = k)
  within_rounding(o$history$statistic[o$history$test == "aux_const"], ex$aux_const_t, 2)
})
