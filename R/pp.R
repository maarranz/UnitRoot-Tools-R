# Phillips-Perron test (spec: Decision D12): internal implementation of the
# Z_tau statistic of Phillips and Perron (1988), in the form of Hamilton (1994,
# eq. 17.6.8), replacing urca::ur.pp (whose correction term uses y_t where the
# published formula uses y_{t-1}).

#' Phillips-Perron test
#'
#' Phillips-Perron \eqn{Z_\tau} test of a unit root. The Dickey-Fuller
#' regression without lagged differences,
#' \deqn{\Delta y_t = a + \gamma t + b y_{t-1} + u_t, \quad t = 2, \dots, T,}
#' (\eqn{\gamma = 0} for `trend = "c"`) is estimated by OLS, and the t ratio
#' of \eqn{b} is corrected for serial correlation in \eqn{u_t}:
#' \deqn{Z_\tau = \sqrt{\hat\gamma_0 / \hat\lambda^2}\, t_b -
#'   \frac{\hat\lambda^2 - \hat\gamma_0}{2 \hat\lambda}\,
#'   \frac{n\, \widehat{se}(\hat b)}{s},}
#' where \eqn{n = T - 1}, \eqn{\hat\gamma_0 = n^{-1} \sum \hat u_t^2},
#' \eqn{s^2} is the OLS residual variance (degrees-of-freedom corrected) and
#' \eqn{\hat\lambda^2} the Bartlett-kernel long-run variance of
#' \eqn{\hat u_t} (Phillips and Perron, 1988; Hamilton, 1994, eq. 17.6.8).
#'
#' Critical values and the p-value are those of the Dickey-Fuller
#' \eqn{\tau} (MacKinnon's response surfaces and MacKinnon, 1996, at the
#' effective sample size, as in [adf()]): \eqn{Z_\tau} has the same limiting
#' distribution. Its finite-sample distribution differs, especially with
#' strong negative moving-average components.
#'
#' **Long-run variance.** As in [kpss()]: kernel `"bartlett"` (default),
#' `"parzen"` or `"qs"`; bandwidth `"short"` (default) or `"long"`
#' (Schwert's rules applied to \eqn{n = T - 1}), a fixed number of
#' autocovariances, `"andrews"` or `"newey-west"` (automatic).
#'
#' The normalized-bias statistic \eqn{Z_\alpha} is reported in
#' `other_statistics`, without critical values.
#'
#' @inheritParams adf
#' @param trend `"c"` or `"ct"` (or an alias accepted by [match_trend()]).
#' @param level Significance level of the reported decision: 0.01, 0.05 or
#'   0.10.
#' @param kernel `"bartlett"` (default), `"parzen"` or `"qs"`.
#' @param bandwidth `"short"` (default), `"long"`, `"andrews"`,
#'   `"newey-west"`, or a non-negative integer number of autocovariances.
#' @return An object of class `ur_test`.
#' @references
#' Hamilton, J. D. (1994). *Time Series Analysis*. Princeton University
#' Press, Section 17.6.
#'
#' Phillips, P. C. B. and Perron, P. (1988). Testing for a unit root in time
#' series regression. *Biometrika*, 75, 335-346.
#' @examples
#' set.seed(1)
#' pp(cumsum(rnorm(150)), trend = "c")
#' @export
pp <- function(y, trend = "c", kernel = "bartlett", bandwidth = "short",
               level = 0.05, var = NULL) {
  call <- match.call()
  check_level(level)
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  kernel <- check_kernel(kernel)
  n <- length(y) - 1L

  r <- adf_regression(y, trend, 0L, 2L)
  f <- r$fit
  u <- f$residuals
  bw <- resolve_bandwidth(bandwidth, u, kernel, n_rule = n, n_label = "(T-1)")
  if (bw$b - 1 >= n) {
    stop("`bandwidth` is too large for T = ", n + 1L, ".", call. = FALSE)
  }
  v <- long_run_variance(u, kernel, bw$b)
  g0 <- v$gamma0
  lrv <- v$lrv
  s <- sqrt(f$sigma2)
  se_b <- f$se[["y_lag1"]]
  stat <- sqrt(g0 / lrv) * f$t_values[["y_lag1"]] -
    0.5 * (lrv - g0) / sqrt(lrv) * n * se_b / s
  # normalized bias (Phillips and Perron 1988; Hamilton 1994, eq. 17.6.7)
  z_alpha <- n * f$coefficients[["y_lag1"]] -
    0.5 * (n^2 * se_b^2 / s^2) * (lrv - g0)

  coef_table <- cbind(Estimate = f$coefficients, `Std. Error` = f$se,
                      `t value` = f$t_values)
  new_ur_test(
    method = "Phillips-Perron test", test = "PP", statistic_name = "Z_tau",
    statistic = unname(stat), null = "unit root",
    alternative = adf_alternatives[[trend]], trend = trend,
    lags = v$lags_used, bandwidth_rule = bandwidth_text(kernel, bw),
    nobs = r$nobs, sample = unname(r$sample),
    critical_values = tau_critical_values(trend, r$nobs), tail = "lower",
    p_value = tau_p_value(stat, trend, r$nobs),
    calibration = list(
      method = "response_surface", source = mackinnon_source[[trend]],
      note = sprintf("evaluated at N = %d; Z_tau has the Dickey-Fuller tau limit distribution", r$nobs),
      p_value_source = "MacKinnon (1996) finite-sample p-value for tau (urca::punitroot)"),
    level = level,
    other_statistics = c(Z_alpha = unname(z_alpha),
                         t_b = unname(f$t_values[["y_lag1"]]),
                         gamma0 = unname(g0), long_run_variance = unname(lrv),
                         bandwidth = unname(v$b)),
    engine = engine_internal(),
    warnings = if (r$nobs < 25L) "N is below 25." else character(),
    regression = list(coefficients = coef_table, rss = f$rss, df = f$df,
                      sigma2 = f$sigma2),
    call = call)
}
