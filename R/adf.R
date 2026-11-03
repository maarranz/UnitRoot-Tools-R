# Augmented Dickey-Fuller test (spec: "The ADF regression")

adf_alternatives <- c(
  n  = "stationary with zero mean",
  c  = "stationary around a constant level",
  ct = "stationary around a linear trend"
)

#' Augmented Dickey-Fuller test
#'
#' Tests the null hypothesis of a unit root, \eqn{H_0: b = 0}, against
#' stationarity, \eqn{H_1: b < 0}, in the regression
#' \deqn{\Delta y_t = a + \gamma t + b y_{t-1} + \sum_{j=1}^{k} \beta_j
#'   \Delta y_{t-j} + \varepsilon_t,}
#' with \eqn{a = \gamma = 0} imposed for `trend = "n"` and \eqn{\gamma = 0}
#' for `trend = "c"`.
#'
#' The statistic is the t ratio of \eqn{b} (\eqn{\tau_1}, \eqn{\tau_2} or
#' \eqn{\tau_3}). Its critical values come from MacKinnon's response surfaces
#' evaluated at the effective sample size. The joint statistics
#' \eqn{\phi_1} (`"c"`) or \eqn{\phi_2} and \eqn{\phi_3} (`"ct"`) are
#' reported alongside, with Dickey and Fuller (1981) critical values; they
#' are computed from restricted regressions on the same sample.
#'
#' **Lags.** Give `lags` to fix the number of lagged differences; the
#' regression then uses the maximal sample, \eqn{t = k + 2, \dots, T}.
#' Otherwise the order is chosen by [select_lags()] over
#' \eqn{0, \dots, k_{\max}} on the common sample \eqn{t = k_{\max} + 2, \dots,
#' T}, and the test regression is estimated on that same sample
#' (`final_sample = "common"`) or re-estimated on the maximal sample for the
#' selected order (`final_sample = "maximal"`).
#'
#' @inheritParams select_lags
#' @param lags Fixed number of lagged differences. If `NULL` (default), the
#'   order is selected by `criterion`.
#' @param final_sample After lag selection, estimate the test regression on
#'   the `"common"` sample used for selection (default) or on the
#'   `"maximal"` sample for the selected order.
#' @param level Significance level of the reported decision: 0.01, 0.05 or
#'   0.10 (the levels with tabulated critical values).
#'
#' @return An object of class `ur_test` (see [new_ur_test()]).
#' @references
#' Dickey, D. A. and Fuller, W. A. (1981). Likelihood ratio statistics for
#' autoregressive time series with a unit root. *Econometrica*, 49, 1057-1072.
#'
#' MacKinnon, J. G. (2010). Critical values for cointegration tests. Queen's
#' Economics Department Working Paper 1227.
#' @examples
#' set.seed(1)
#' y <- cumsum(0.2 + rnorm(150))
#' adf(y, trend = "ct")
#' adf(y, trend = "ct", lags = 2)
#' summary(adf(y, trend = "c"))
#' @export
adf <- function(y, trend = "c", lags = NULL, max_lags = NULL,
                criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                final_sample = c("common", "maximal"),
                level = 0.05, var = NULL) {
  call <- match.call()
  criterion <- match.arg(criterion)
  final_sample <- match.arg(final_sample)
  check_level(level)
  y <- check_series(y, var = var)
  trend <- match_trend(trend)

  if (!is.null(lags)) {
    if (!is.null(max_lags)) {
      stop("Give either `lags` (fixed order) or `max_lags` (selection), not both.",
           call. = FALSE)
    }
    k <- check_order(lags, "lags")
    first <- k + 2L
    selection <- NULL
  } else {
    sel <- lag_table(y, trend, max_lags, criterion)
    k <- sel$lags
    first <- if (final_sample == "common") sel$kmax + 2L else k + 2L
    selection <- list(method = criterion, kmax = sel$kmax,
                      detrending = sel$detrending, table = sel$table,
                      final_sample = final_sample)
  }

  r <- adf_regression(y, trend, k, first)
  adf_result(r, selection, level, call)
}

# Builds the ur_test for an ADF regression (used by adf() and the strategies).
adf_result <- function(r, selection, level, call) {
  trend <- r$trend
  warn <- character()
  if (r$nobs < 25L) {
    warn <- c(warn, sprintf(
      "N = %d is below 25, the smallest sample size in the critical-value tables.",
      r$nobs))
  }

  phi_cv <- lapply(names(r$phi), function(nm) phi_critical_values(nm, r$nobs))
  names(phi_cv) <- names(r$phi)
  if (any(vapply(phi_cv, `[[`, logical(1), "extrapolated"))) {
    warn <- c(warn, "phi critical values: T = 25 row used for N < 25 (extrapolation).")
  }

  fit <- r$fit
  coef_table <- cbind(Estimate = fit$coefficients, `Std. Error` = fit$se,
                      `t value` = fit$t_values)

  new_ur_test(
    method = "Augmented Dickey-Fuller test",
    test = "ADF",
    statistic_name = names(r$tau),
    statistic = unname(r$tau),
    null = "unit root",
    alternative = adf_alternatives[[trend]],
    trend = trend,
    lags = r$k,
    lag_selection = selection,
    nobs = r$nobs,
    sample = unname(r$sample),
    critical_values = tau_critical_values(trend, r$nobs),
    tail = "lower",
    calibration = list(
      method = "response_surface",
      source = mackinnon_source[[trend]],
      note = sprintf("evaluated at N = %d", r$nobs)
    ),
    level = level,
    other_statistics = r$phi,
    other_critical_values = lapply(phi_cv, function(p) {
      structure(p$values, source = sprintf(
        "Dickey and Fuller (1981), table row T = %s", format(p$size)))
    }),
    engine = engine_internal(),
    warnings = warn,
    regression = list(coefficients = coef_table, rss = fit$rss,
                      df = fit$df, sigma2 = fit$sigma2),
    call = call
  )
}
