# Adapters to urca (spec: "Known backend limitations")
#
# Each function translates arguments and conventions, calls urca, and returns
# the common `ur_test` result. Only options urca actually supports are
# exposed. The `engine` field records the backend so that results can be
# compared when these adapters are replaced by internal implementations.

engine_urca <- function() {
  paste("urca", as.character(utils::packageVersion("urca")))
}

#' DF-GLS test
#'
#' Elliott-Rothenberg-Stock DF-GLS test of a unit root: the series is
#' detrended by GLS (local-to-unity, \eqn{\bar c = -7} or \eqn{-13.5}) and a
#' Dickey-Fuller regression without deterministic terms is run on the
#' detrended series. Computed by [urca::ur.ers()].
#'
#' **Lags.** Give `lags` to fix the number of lagged differences. Otherwise
#' the order is selected as recommended by Perron and Qu (2007): by
#' [select_lags()] (MAIC by default) on the series detrended by **OLS**,
#' exactly as for [adf()]; the DF-GLS statistic is then computed on the
#' GLS-detrended series with that number of lags, on the maximal sample.
#'
#' **Critical values and p-value.** The DF-GLS statistic has a heavier
#' left tail in finite samples than the Dickey-Fuller distributions often
#' used for it, and its quantiles also depend on the number of lagged
#' differences. Critical values and the p-value therefore come from this
#' package's own response surfaces, estimated from a Monte Carlo simulation
#' of the statistic under the null for \eqn{T = 25, \dots, 1000} and
#' \eqn{k = 0, \dots, 21} (see the specification, Decision D11). A warning
#' is given when \eqn{(T, k)} lies outside the simulated range.
#'
#' The statistic itself is computed by [urca::ur.ers()] (checked against an
#' independent implementation).
#'
#' @inheritParams pp
#' @inheritParams select_lags
#' @param lags Number of lagged differences in the DF-GLS regression. If
#'   `NULL` (default), selected by `criterion` as described above.
#' @return An object of class `ur_test`.
#' @references
#' Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996). Efficient tests for
#' an autoregressive unit root. *Econometrica*, 64, 813-836.
#'
#' Ng, S. and Perron, P. (2001). Lag length selection and the construction of
#' unit root tests with good size and power. *Econometrica*, 69, 1519-1554.
#'
#' Perron, P. and Qu, Z. (2007). A simple modification to improve the finite
#' sample properties of Ng and Perron's unit root tests. *Economics Letters*,
#' 94, 12-19.
#' @examples
#' set.seed(1)
#' y <- cumsum(rnorm(150))
#' dfgls(y, trend = "c")
#' dfgls(y, trend = "c", lags = 2)
#' @export
dfgls <- function(y, trend = "c", lags = NULL, max_lags = NULL,
                  criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                  level = 0.05, var = NULL) {
  call <- match.call()
  criterion <- match.arg(criterion)
  check_level(level)
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  notes <- character()
  if (!is.null(lags)) {
    if (!is.null(max_lags)) {
      stop("Give either `lags` (fixed order) or `max_lags` (selection), not both.",
           call. = FALSE)
    }
    k <- check_order(lags, "lags")
    selection <- NULL
  } else {
    sel <- lag_table(y, trend, max_lags, criterion)
    k <- sel$lags
    notes <- sel$notes
    selection <- list(method = criterion, kmax = sel$kmax,
                      detrending = "OLS (Perron and Qu, 2007); statistic on GLS-detrended series",
                      table = sel$table, final_sample = "maximal", notes = sel$notes)
  }
  if (length(y) - k - 1L < k + 2L) {
    stop("`lags` = ", k, " is too large for T = ", length(y), ".", call. = FALSE)
  }
  u <- urca::ur.ers(y, type = "DF-GLS",
                    model = if (trend == "c") "constant" else "trend",
                    lag.max = k)
  n <- length(y)
  nobs <- n - k - 1L
  cv <- dfgls_critical_values(trend, n, k)
  calibration <- list(
    method = "response_surface",
    source = "unitroottools DF-GLS response surfaces (own simulation)",
    note = sprintf("finite-sample quantiles for T = %d and k = %d lagged differences", n, k))
  design_warnings <- dfgls_design_warnings(trend, n, k)
  new_ur_test(
    method = "DF-GLS test", test = "DF-GLS",
    statistic_name = "tau_GLS", statistic = unname(u@teststat),
    null = "unit root", alternative = adf_alternatives[[trend]], trend = trend,
    lags = k, lag_selection = selection, nobs = nobs, sample = c(k + 2L, n),
    critical_values = cv, tail = "lower",
    p_value = dfgls_p_value(unname(u@teststat), trend, n, k),
    calibration = calibration, level = level, engine = engine_urca(),
    warnings = c(notes, design_warnings),
    call = call)
}
