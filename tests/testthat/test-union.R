# union_test(): Harvey, Leybourne and Taylor (2012) union of rejections with
# own finite-sample calibration (D22).

test_that("the components are adf() and dfgls() with the common lag order", {
  set.seed(31)
  y <- cumsum(rnorm(180))
  for (k in c(0L, 3L)) {
    u <- union_test(y, lags = k)
    expect_equal(u$other_statistics[["ADF_c"]], adf(y, trend = "c", lags = k)$statistic, tolerance = 1e-12)
    expect_equal(u$other_statistics[["ADF_ct"]], adf(y, trend = "ct", lags = k)$statistic, tolerance = 1e-12)
    expect_equal(u$other_statistics[["GLS_c"]], dfgls(y, trend = "c", lags = k)$statistic, tolerance = 1e-12)
    expect_equal(u$other_statistics[["GLS_ct"]], dfgls(y, trend = "ct", lags = k)$statistic, tolerance = 1e-12)
    # each component carries the critical values and p-value of its own test
    expect_equal(as.numeric(u$other_critical_values$GLS_ct), unname(dfgls(y, trend = "ct", lags = k)$critical_values))
    expect_equal(attr(u$other_critical_values$ADF_c, "p_value"), adf(y, trend = "c", lags = k)$p_value)
    expect_identical(u$nobs, length(y) - k - 1L)
    expect_identical(u$sample, c(k + 2L, length(y)))
  }
})

test_that("U is the smallest normalised statistic; decision and p-value from the union surface", {
  set.seed(32)
  y <- cumsum(rnorm(150))
  u <- union_test(y, lags = 2)
  cv5 <- union_cv5(150, 2)
  z <- u$other_statistics / abs(cv5)
  expect_equal(u$statistic, min(z))
  expect_match(u$calibration$note, names(z)[which.min(z)], fixed = TRUE)
  expect_equal(u$critical_values, rs_critical_values(union_rs, 150, 2))
  expect_equal(u$p_value, rs_p_value(union_rs, u$statistic, 150, 2))
  expect_identical(u$reject, u$statistic < u$critical_values[["0.05"]])
  # the union critical value is beyond -1 (each component's own 5% point)
  expect_lt(u$critical_values[["0.05"]], -1)
  expect_gt(u$critical_values[["0.05"]], -1.5)
})

test_that("normalisation: each component's 5% critical value maps to -1", {
  for (T in c(40L, 200L)) for (k in c(0L, 4L)) {
    cv5 <- union_cv5(T, k)
    expect_equal(unname(cv5[["ADF_ct"]]),
                 unname(tau_calibration("ct", 0, T - k - 1L, k, levels = 0.05)$critical_values))
    expect_equal(unname(cv5[["GLS_c"]]), unname(dfgls_critical_values("c", T, k)[["0.05"]]))
  }
})

test_that("stationary series are rejected; lag selection, input checks", {
  set.seed(33)
  e <- rnorm(250)
  u <- union_test(e)
  expect_true(u$reject)
  expect_identical(u$lag_selection$method, "MAIC")
  expect_identical(union_test(e, criterion = "bic")$lag_selection$method, "BIC")
  expect_error(union_test(rnorm(30), lags = 13), "too large")
  expect_error(union_test(e, level = 0.02), "level")
})

test_that("residual checks, tidy and printing work on union results", {
  set.seed(34)
  y <- cumsum(rnorm(120))
  u <- union_test(y, lags = 1)
  expect_s3_class(residual_check(u), "ur_residual_check")
  td <- tidy(u)
  expect_identical(td$statistic_name, c("U", "ADF_c", "ADF_ct", "GLS_c", "GLS_ct"))
  expect_identical(glance(u)$test, "Union")
  expect_output(print(u), "Union-of-rejections")
  expect_output(print(u), "constant \\(2 tests\\) and constant and trend \\(2 tests\\)")
  expect_output(print(summary(u)), "GLS_ct")
})

test_that("the shared union surface (JSON) is the package's (D22)", {
  skip_if_not_installed("jsonlite")
  j <- jsonlite::read_json(system.file("fixtures", "union_response_surface.json",
                                       package = "unitroottools"))
  jt <- j$trends$union
  expect_identical(unlist(jt$intercept), union_rs$intercept)
  for (a in seq_along(union_rs$ks)) {
    expect_identical(unlist(jt$slopes[[a]][[5L]]), unname(union_rs$slopes[a, 5L, ]))
  }
  for (row in j$reference_values$rows) {
    expect_equal(unname(rs_critical_values(union_rs, row$T, row$k)), c(row$cv01, row$cv05, row$cv10),
                 tolerance = 1e-12)
    expect_equal(rs_p_value(union_rs, row$stat, row$T, row$k), row$p, tolerance = 1e-12)
    expect_equal(unname(union_cv5(row$T, row$k)),
                 c(row$cv5_ADF_c, row$cv5_ADF_ct, row$cv5_GLS_c, row$cv5_GLS_ct), tolerance = 1e-12)
  }
})
