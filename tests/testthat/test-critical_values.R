test_that("tau critical values reproduce MacKinnon's response surfaces", {
  # Reference values from statsmodels 0.15 (adfvalues.mackinnoncrit, nobs = 100)
  expect_equal(unname(tau_critical_values("c", 100)),
               c(-3.49750103, -2.89090644, -2.5824349), tolerance = 1e-8)
  expect_equal(unname(tau_critical_values("ct", 100)),
               c(-4.05227796, -3.45534297, -3.15332088), tolerance = 1e-8)
  expect_equal(unname(tau_critical_values("n", 100)),
               c(-2.5884607, -1.94399128, -1.61441004), tolerance = 1e-8)
  # Asymptotic limits are the leading coefficients
  expect_equal(unname(tau_critical_values("c", 1e12)),
               c(-3.43035, -2.86154, -2.56677), tolerance = 1e-8)
  expect_named(tau_critical_values("c", 50), c("0.01", "0.05", "0.1"))
})

test_that("phi table rows: largest tabulated size not exceeding N, capped at 500", {
  expect_identical(phi_table_row(49)$size, 25)
  expect_identical(phi_table_row(50)$size, 50)
  expect_identical(phi_table_row(249)$size, 100)
  expect_identical(phi_table_row(250)$size, 250)
  expect_identical(phi_table_row(10000)$size, 500)
  r <- phi_table_row(20)
  expect_identical(r$size, 25)
  expect_true(r$extrapolated)
})

test_that("phi critical values are the Dickey-Fuller (1981) table entries", {
  expect_equal(unname(phi_critical_values("phi1", 100)$values), c(6.70, 4.71, 3.86))
  expect_equal(unname(phi_critical_values("phi2", 100)$values), c(6.50, 4.88, 4.16))
  expect_equal(unname(phi_critical_values("phi3", 100)$values), c(8.73, 6.49, 5.47))
  # Published T = 250 row of phi3 (urca 1.3-4 repeats the T = 100 values here)
  expect_equal(unname(phi_critical_values("phi3", 300)$values), c(8.43, 6.34, 5.39))
})
