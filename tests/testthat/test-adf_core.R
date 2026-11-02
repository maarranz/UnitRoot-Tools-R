# The ADF core is checked against independent lm() fits on explicitly built
# data. Tolerance follows the specification's numerical-reproducibility
# standard (relative 1e-10).

tol <- 1e-10

sim_rw <- function(n = 120, seed = 1, drift = 0.1) {
  set.seed(seed)
  cumsum(drift + stats::arima.sim(list(ar = 0.4), n = n))
}

# Independent reference: build the regression with embed() and fit with lm()
lm_reference <- function(y, trend, k, first) {
  n <- length(y)
  t <- first:n
  dy <- diff(y)                                   # dy[i] = y[i+1] - y[i]
  df <- data.frame(dy = dy[t - 1], ylag = y[t - 1], tt = t)
  for (j in seq_len(k)) df[[paste0("l", j)]] <- dy[t - 1 - j]
  lags <- if (k > 0) paste0("l", seq_len(k)) else character()
  rhs <- c(switch(trend, n = "0", c = character(), ct = "tt"), "ylag", lags)
  f <- stats::as.formula(paste("dy ~", paste(rhs, collapse = " + ")))
  list(fit = stats::lm(f, data = df), data = df, lags = lags)
}

rss_lm <- function(f, data) sum(stats::resid(stats::lm(f, data = data))^2)

test_that("design matrix uses the original time index and the requested sample", {
  y <- sim_rw()
  d <- adf_design(y, "ct", k = 3, first = 10)
  expect_identical(d$t, 10:120)
  expect_identical(colnames(d$X), c("const", "trend", "y_lag1",
                                    "dy_lag1", "dy_lag2", "dy_lag3"))
  expect_equal(d$X[, "trend"], as.numeric(10:120))
  expect_equal(d$response[1], y[10] - y[9])
  expect_equal(unname(d$X[1, "y_lag1"]), y[9])
  expect_equal(unname(d$X[1, "dy_lag3"]), y[7] - y[6])
  expect_identical(colnames(adf_design(y, "n", 0)$X), "y_lag1")
  expect_identical(adf_design(y, "c", 2)$t[1], 4L)       # maximal sample: k + 2
})

test_that("invalid samples are rejected", {
  y <- sim_rw()
  expect_error(adf_design(y, "c", k = 3, first = 4), "t >= k \\+ 2 = 5")
  expect_error(adf_design(y, "c", k = -1), "non-negative")
  expect_error(adf_design(y, "c", k = 0, first = 121), "empty")
})

test_that("coefficients, standard errors, t ratios and RSS match lm()", {
  y <- sim_rw()
  for (trend in c("n", "c", "ct")) for (k in c(0, 2, 5)) {
    first <- 9                                       # common sample, kmax = 7
    ref <- lm_reference(y, trend, k, first)
    s <- summary(ref$fit)$coefficients
    r <- adf_regression(y, trend, k, first)
    expect_equal(unname(r$fit$coefficients), unname(s[, 1]), tolerance = tol)
    expect_equal(unname(r$fit$se), unname(s[, 2]), tolerance = tol)
    expect_equal(unname(r$fit$t_values), unname(s[, 3]), tolerance = tol)
    expect_equal(r$fit$rss, sum(stats::resid(ref$fit)^2), tolerance = tol)
    expect_equal(unname(r$tau), s["ylag", 3], tolerance = tol)
    expect_identical(r$nobs, 112L)
    expect_identical(unname(r$sample), c(9L, 120L))
  }
})

test_that("tau carries the name of its deterministic specification", {
  y <- sim_rw()
  expect_named(adf_regression(y, "n", 1)$tau, "tau1")
  expect_named(adf_regression(y, "c", 1)$tau, "tau2")
  expect_named(adf_regression(y, "ct", 1)$tau, "tau3")
})

test_that("phi statistics match F tests of nested lm() fits on the same sample", {
  y <- sim_rw(seed = 7)
  for (k in c(0, 3)) {
    first <- 6
    ref <- lm_reference(y, "ct", k, first)
    dat <- ref$data
    lags <- if (length(ref$lags)) paste(ref$lags, collapse = " + ") else NULL
    f_u_ct <- stats::as.formula(paste("dy ~", paste(c("tt", "ylag", lags), collapse = " + ")))
    f_phi3 <- stats::as.formula(paste("dy ~", paste(c("1", lags), collapse = " + ")))
    f_phi2 <- stats::as.formula(paste("dy ~", paste(c("0", lags), collapse = " + ")))
    ct <- adf_regression(y, "ct", k, first)
    expect_equal(ct$phi[["phi3"]],
                 stats::anova(stats::lm(f_phi3, dat), stats::lm(f_u_ct, dat))$F[2],
                 tolerance = tol)
    expect_equal(ct$phi[["phi2"]],
                 stats::anova(stats::lm(f_phi2, dat), stats::lm(f_u_ct, dat))$F[2],
                 tolerance = tol)
    expect_named(ct$phi, c("phi2", "phi3"))

    f_u_c <- stats::as.formula(paste("dy ~", paste(c("ylag", lags), collapse = " + ")))
    cc <- adf_regression(y, "c", k, first)
    expect_equal(cc$phi[["phi1"]],
                 stats::anova(stats::lm(f_phi2, dat), stats::lm(f_u_c, dat))$F[2],
                 tolerance = tol)
    expect_named(cc$phi, "phi1")
    expect_length(adf_regression(y, "n", k, first)$phi, 0)
  }
})

test_that("a restricted model with no regressors is handled (phi1 with k = 0)", {
  y <- sim_rw(seed = 3)
  r <- adf_regression(y, "c", 0)
  dy <- diff(y)
  rss_r <- sum(dy^2)
  expect_equal(r$phi[["phi1"]],
               ((rss_r - r$fit$rss) / 2) / (r$fit$rss / (r$nobs - 2)),
               tolerance = tol)
})

test_that("degenerate regressions give informative errors", {
  expect_error(adf_regression(rep(1, 30), "c", 0), "collinear")
  expect_error(adf_regression(cumsum(rnorm(6)), "ct", 3), "Not enough observations")
})
