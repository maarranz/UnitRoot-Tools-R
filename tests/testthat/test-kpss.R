# KPSS robustness (audit of 0.2.0, findings A1 and N2)

test_that("large statistics return a result with a bounded p-value (A1)", {
  r <- kpss(seq_len(1000), bandwidth = 0)
  expect_s3_class(r, "ur_test")
  expect_true(r$reject)
  expect_equal(r$p_value, 0.001)
  expect_match(r$calibration$note, "below 0.001")
  expect_match(capture.output(print(r)), "p-value: +< 0.001", all = FALSE)
  # p-values: monotone over the whole range, bounded, right at the critical values
  for (tr in c("c", "ct")) {
    s <- c(seq(0.001, 3, length.out = 120), 10, 100)
    p <- vapply(s, kpss_p_value, numeric(1), trend = tr)
    expect_true(all(diff(p) <= 1e-12))
    expect_true(all(p >= 0.001 & p <= 0.999))
    cv <- kpss_critical_values(tr)
    expect_equal(vapply(cv, kpss_p_value, numeric(1), trend = tr),
                 c(`0.1` = 0.10, `0.05` = 0.05, `0.025` = 0.025, `0.01` = 0.01), tolerance = 1e-6)
  }
})

test_that("constant series and exact trends stop with a clear error (N2)", {
  expect_error(kpss(rep(1, 100)), "the series is constant")
  expect_error(kpss(rep(0, 100)), "the series is constant")
  expect_error(kpss(3 + 0.5 * seq_len(100), trend = "ct"), "exact linear trend")
  expect_s3_class(kpss(3 + 0.5 * seq_len(100) + rnorm(100), trend = "ct"), "ur_test")
})
