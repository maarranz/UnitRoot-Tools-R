test_that("simulated quantiles are close to MacKinnon's critical values", {
  s <- simulate_df(T = 100, R = 4000, trend = "c", seed = 123)
  expect_s3_class(s, "df_simulation")
  expect_length(s$tau, 4000)
  q <- s$quantiles
  # Monte Carlo error of a 5% quantile with R = 4000 is about 0.03
  expect_lt(abs(q$simulated[q$level == 0.05] - q$mackinnon[q$level == 0.05]), 0.1)
  expect_lt(q$simulated[q$level == 0.05], stats::qnorm(0.05))
})

test_that("simulations are reproducible and restore the RNG state", {
  set.seed(99); before <- stats::runif(1)
  set.seed(99)
  a <- simulate_df(T = 50, R = 200, seed = 5)
  after <- stats::runif(1)
  expect_identical(after, before)
  b <- simulate_df(T = 50, R = 200, seed = 5)
  expect_identical(a$tau, b$tau)
})

test_that("tau equals the ADF statistic on the same series", {
  set.seed(7)
  y <- cumsum(stats::rnorm(60))
  set.seed(7)
  s <- simulate_df(T = 60, R = 100, trend = "ct")
  expect_equal(s$tau[1], adf(y, trend = "ct", lags = 0)$statistic, tolerance = 1e-12)
})

test_that("arguments are checked; print and plot work", {
  expect_error(simulate_df(T = 5), "at least 10")
  expect_error(simulate_df(R = 10), "at least 100")
  s <- simulate_df(T = 50, R = 200, seed = 1)
  out <- capture.output(print(s))
  expect_match(out, "Simulated Dickey-Fuller distribution of tau2", all = FALSE)
  expect_match(out, "MacKinnon", all = FALSE)
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(s))
})
