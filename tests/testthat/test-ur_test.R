# A valid ADF-type result; individual tests override fields.
make_result <- function(...) {
  args <- list(
    method = "Augmented Dickey-Fuller test",
    test = "ADF",
    statistic_name = "tau3",
    statistic = -3.80,
    null = "unit root",
    alternative = "stationary around a linear trend",
    trend = "ct",
    lags = 4L,
    lag_selection = list(method = "MAIC", kmax = 12,
                         detrending = "OLS",
                         table = data.frame(k = 0:12, MAIC = seq(-1, -2.2, by = -0.1))),
    nobs = 187L,
    sample = c(14L, 200L),
    critical_values = c("0.01" = -4.00, "0.05" = -3.43, "0.1" = -3.14),
    tail = "lower",
    calibration = list(method = "response_surface", source = "MacKinnon (2010)"),
    other_statistics = c(phi2 = 5.1, phi3 = 7.4),
    engine = "internal 0.0.0.9000"
  )
  args[names(list(...))] <- list(...)
  do.call(new_ur_test, args)
}

test_that("a valid result is constructed with a decision", {
  r <- make_result()
  expect_s3_class(r, "ur_test")
  expect_true(r$reject)                      # -3.80 < -3.43
  expect_identical(r$trend, "ct")
  expect_identical(r$lags, 4L)
})

test_that("lower-tail decisions compare against the critical value", {
  expect_false(make_result(statistic = -3.00)$reject)
  expect_true(make_result(statistic = -3.20, level = 0.1)$reject)
  expect_false(make_result(statistic = -3.80, level = 0.01)$reject)
})

test_that("upper-tail decisions (KPSS-type) compare in the other direction", {
  kpss <- make_result(
    test = "KPSS", statistic_name = "eta_tau", statistic = 0.20,
    null = "stationary around a linear trend", alternative = "unit root",
    lags = 4L, lag_selection = NULL, bandwidth_rule = "Schwert short",
    critical_values = c("0.01" = 0.216, "0.05" = 0.146, "0.1" = 0.119),
    tail = "upper",
    calibration = list(method = "table", source = "KPSS (1992), Table 1"),
    other_statistics = numeric()
  )
  expect_true(kpss$reject)
  expect_false(make_result(
    statistic = 0.10, tail = "upper", lag_selection = NULL,
    critical_values = c("0.01" = 0.216, "0.05" = 0.146, "0.1" = 0.119)
  )$reject)
})

test_that("a level without a critical value uses the p-value", {
  r <- make_result(level = 0.025, p_value = 0.02)
  expect_true(r$reject)
  expect_error(make_result(level = 0.025), "decision cannot be made")
})

test_that("aliases are accepted and stored canonically", {
  expect_identical(make_result(trend = "trend")$trend, "ct")
})

test_that("invalid contents are rejected", {
  expect_error(make_result(statistic = NA_real_), "finite number")
  expect_error(make_result(null = ""), "non-empty string")
  expect_error(make_result(nobs = 100L), "length of `sample`")
  expect_error(make_result(sample = c(200L, 14L)), "first <= last")
  expect_error(make_result(lags = -1L), "non-negative")
  expect_error(make_result(lags = 13L), "kmax")
  expect_error(make_result(critical_values = c(-4.00, -3.43, -3.14)),
               "named numeric")
  expect_error(make_result(critical_values = c("0.01" = -4, "0.05" = -3.4)),
               "0.01, 0.05 and 0.10")
  expect_error(make_result(critical_values = c("0.01" = -3.14, "0.05" = -3.43, "0.1" = -4)),
               "lower-tail")
  expect_error(make_result(p_value = 1.5), "p_value")
  expect_error(make_result(calibration = list(method = "interpolated", source = "x")),
               "calibration")
  expect_error(make_result(other_statistics = c(5.1, 7.4)), "named numeric")
  expect_error(make_result(tail = "two"), "should be one of")
})

test_that("print shows hypotheses, sample, critical values and decision", {
  out <- capture.output(print(make_result()))
  expect_match(out, "Augmented Dickey-Fuller test \\(tau3\\)", all = FALSE)
  expect_match(out, "H0: unit root", all = FALSE)
  expect_match(out, "Deterministics: constant and linear trend", all = FALSE)
  expect_match(out, "Lags: 4 \\(MAIC, kmax = 12\\)", all = FALSE)
  expect_match(out, "N = 187 \\(t = 14, ..., 200\\)", all = FALSE)
  expect_match(out, "1%: -4.000   5%: -3.430   10%: -3.140", all = FALSE)
  expect_match(out, "Decision at 5%: reject H0 \\(unit root\\)", all = FALSE)
  expect_false(any(grepl("p-value", out)))   # table/response surface without p-value
})

test_that("print reports fixed lags, bandwidths and warnings", {
  out <- capture.output(print(make_result(lag_selection = NULL,
                                          warnings = "N < 25: critical values extrapolated")))
  expect_match(out, "Lags: 4 \\(fixed\\)", all = FALSE)
  expect_match(out, "Warning: N < 25", all = FALSE)
})

test_that("summary adds calibration, other statistics, lag table and engine", {
  out <- capture.output(print(summary(make_result())))
  expect_match(out, "Source: MacKinnon \\(2010\\)", all = FALSE)
  expect_match(out, "not provided by this calibration", all = FALSE)
  expect_match(out, "phi3", all = FALSE)
  expect_match(out, "Lag selection: MAIC, kmax = 12, detrending: OLS", all = FALSE)
  expect_match(out, "Engine: internal 0.0.0.9000", all = FALSE)
})

test_that("M3: the validator enforces the whole contract", {
  r <- make_result()
  bad <- r; bad$reject <- !bad$reject
  expect_error(validate_ur_test(bad), "stored decision")
  bad <- r; bad$reject <- NA
  expect_error(validate_ur_test(bad), "`reject` must be TRUE or FALSE")
  bad <- r; bad$trend <- "drift"
  expect_error(validate_ur_test(bad), "canonical codes")
  bad <- r; bad$tail <- "both"
  expect_error(validate_ur_test(bad), "`tail` must be")
  bad <- r; bad$lags <- 4.5
  expect_error(validate_ur_test(bad), "`lags` must be a non-negative integer")
  bad <- r; bad$nobs <- 187
  expect_error(validate_ur_test(bad), "`nobs` must be a positive integer")
  bad <- r; bad$other_statistics[["phi3"]] <- NaN
  expect_error(validate_ur_test(bad), "finite values")
  bad <- r; bad$call <- "adf(y)"
  expect_error(validate_ur_test(bad), "`call` must be NULL or a call")
  bad <- r; class(bad) <- "list"
  expect_error(validate_ur_test(bad), "not an object of class")
  expect_invisible(validate_ur_test(r))
})

test_that("M3: the constructor never truncates whole-number fields", {
  expect_error(make_result(lags = 4.5), "`lags` must be a whole number or NA")
  expect_error(make_result(nobs = 187.2), "`nobs` must be a whole number")
  expect_error(make_result(sample = c(14.5, 200)), "`sample` must be 2 whole numbers")
  expect_identical(make_result(lags = 4)$lags, 4L)
})

test_that("M3: accompanying critical values get the same checks as the main ones", {
  ocv <- list(phi3 = c("0.01" = 8.73, "0.05" = 6.49, "0.1" = 5.47))
  expect_s3_class(make_result(other_critical_values = ocv), "ur_test")
  ocv_bad <- list(phi3 = c("0.01" = 5.47, "0.05" = 6.49, "0.1" = 8.73))
  expect_error(make_result(other_critical_values = ocv_bad), "upper-tail")
  ocv_bad <- list(phi3 = c("0.05" = 6.49))
  expect_error(make_result(other_critical_values = ocv_bad), "0.01, 0.05 and 0.10")
})

test_that("print shows the source of the critical values", {
  out <- capture.output(print(make_result()))
  expect_match(out, "Source:\\s+MacKinnon \\(2010\\)", all = FALSE)
})
