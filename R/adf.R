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
#' reported alongside; they are computed from restricted regressions on the
#' same sample. Their critical values and p-values come from the package's
#' own finite-sample response surfaces in the sample size and the number of
#' lags (they agree with Dickey and Fuller's 1981 tables without lags).
#'
#' **Lags.** Give `lags` to fix the number of lagged differences; the
#' regression then uses the maximal sample, \eqn{t = k + 2, \dots, T}.
#' Otherwise the order is chosen by [select_lags()] over
#' \eqn{0, \dots, k_{\max}} on the common sample \eqn{t = k_{\max} + 2, \dots,
#' T}, and the test regression is re-estimated on the maximal sample for the
#' selected order, \eqn{t = k^* + 2, \dots, T} (`final_sample = "maximal"`,
#' the default), or kept on the common sample used for selection
#' (`final_sample = "common"`), so that it is exactly the selected row of the
#' selection table.
#'
#' **p-value.** MacKinnon's (1996) finite-sample p-value at the effective
#' sample size, computed by [urca::punitroot()]. The decision is taken from
#' the critical value; the two sources agree to four decimals at the
#' tabulated levels.
#'
#' @inheritParams select_lags
#' @param lags Fixed number of lagged differences. If `NULL` (default), the
#'   order is selected by `criterion`.
#' @param final_sample After lag selection, re-estimate the test regression
#'   on the `"maximal"` sample for the selected order (default) or keep the
#'   `"common"` sample used for selection.
#' @param level Significance level of the reported decision: 0.01, 0.05 or
#'   0.10 (the levels with tabulated critical values).
#'
#' @return An object of class `ur_test` (see [new_ur_test()]).
#' @references
#' Dickey, D. A. and Fuller, W. A. (1981). Likelihood ratio statistics for
#' autoregressive time series with a unit root. *Econometrica*, 49, 1057-1072.
#'
#' MacKinnon, J. G. (1996). Numerical distribution functions for unit root
#' and cointegration tests. *Journal of Applied Econometrics*, 11, 601-618.
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
                final_sample = c("maximal", "common"),
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
                      final_sample = final_sample, notes = sel$notes)
  }

  r <- adf_regression(y, trend, k, first)
  adf_result(r, selection, level, call)
}

# Builds the ur_test for an ADF regression (used by adf() and the strategies).
adf_result <- function(r, selection, level, call) {
  trend <- r$trend
  warn <- if (is.null(selection$notes)) character() else selection$notes
  if (r$nobs < 25L) {
    warn <- c(warn, sprintf(
      "N = %d is below 25, the smallest sample size in the critical-value tables.",
      r$nobs))
  }

  # phi statistics: own response surfaces (D16), indexed by the length of a
  # series whose maximal sample is the effective sample, T = N + k + 1
  T_eq <- r$nobs + r$k + 1L
  phi_cv <- lapply(names(r$phi), function(nm) {
    phi_calibration(nm, r$phi[[nm]], T_eq, r$k)
  })
  names(phi_cv) <- names(r$phi)
  if (length(phi_cv)) {
    w <- unique(unlist(lapply(names(r$phi), function(nm) {
      rs_design_warnings(phi_rs[[nm]], T_eq, r$k)
    })))
    if (length(w)) warn <- c(warn, paste("phi statistics:", w))
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
    p_value = tau_p_value(r$tau, trend, r$nobs),
    calibration = list(
      method = "response_surface",
      source = mackinnon_source[[trend]],
      note = sprintf("evaluated at N = %d", r$nobs),
      p_value_source = "MacKinnon (1996) finite-sample p-value (urca::punitroot)"
    ),
    level = level,
    other_statistics = r$phi,
    other_critical_values = phi_cv,
    engine = engine_internal(),
    warnings = warn,
    regression = list(coefficients = coef_table, rss = fit$rss,
                      df = fit$df, sigma2 = fit$sigma2),
    call = call
  )
}
