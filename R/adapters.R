# Adapters to urca (spec: "Known backend limitations")
#
# Each function translates arguments and conventions, calls urca, and returns
# the common `ur_test` result. Only options urca actually supports are
# exposed. The `engine` field records the backend so that results can be
# compared when these adapters are replaced by internal implementations.

engine_urca <- function() {
  paste("urca", as.character(utils::packageVersion("urca")))
}

bandwidth_rules <- c(
  short = "short: trunc(4 (n/100)^(1/4))",
  long  = "long: trunc(12 (n/100)^(1/4))"
)

# "short", "long" or a non-negative integer -> arguments for urca
bandwidth_args <- function(bandwidth) {
  if (is.character(bandwidth) && length(bandwidth) == 1L &&
      bandwidth %in% names(bandwidth_rules)) {
    return(list(lags = bandwidth, use.lag = NULL,
                rule = bandwidth_rules[[bandwidth]]))
  }
  if (is.numeric(bandwidth) && !is.logical(bandwidth)) {
    return(list(lags = "short", use.lag = check_order(bandwidth, "bandwidth"),
                rule = "fixed by the user"))
  }
  stop("`bandwidth` must be \"short\", \"long\" or a non-negative integer. ",
       "Automatic (Newey-West, Andrews) bandwidths and kernels other than ",
       "Bartlett are not available in this version.", call. = FALSE)
}

#' Phillips-Perron test
#'
#' Phillips-Perron \eqn{Z_\tau} test of a unit root, with the long-run
#' variance estimated by a Bartlett kernel. Computed by [urca::ur.pp()].
#'
#' Critical values are MacKinnon's response surfaces for \eqn{\tau} evaluated
#' at the effective sample size (as in [adf()]): \eqn{Z_\tau} has the same
#' limiting distribution as the Dickey-Fuller \eqn{\tau}.
#'
#' Limitations of this version (backend: urca): only \eqn{Z_\tau} (urca does
#' not provide critical values for \eqn{Z_\alpha}); constant or constant and
#' trend; Bartlett kernel with bandwidth `"short"`, `"long"` (Schwert's rules,
#' applied to \eqn{n = T - 1}) or fixed. urca's correction term uses sums
#' over \eqn{y_t} where Phillips and Perron (1988) use \eqn{y_{t-1}}; the
#' difference vanishes as \eqn{T} grows but reaches about 0.1 at
#' \eqn{T = 30}. An internal implementation is planned.
#'
#' @inheritParams adf
#' @param trend `"c"` or `"ct"` (or an alias accepted by [match_trend()]).
#' @param level Significance level of the reported decision: 0.01, 0.05 or
#'   0.10.
#' @param bandwidth `"short"` (default), `"long"`, or a non-negative integer.
#' @return An object of class `ur_test`.
#' @references
#' Phillips, P. C. B. and Perron, P. (1988). Testing for a unit root in time
#' series regression. *Biometrika*, 75, 335-346.
#' @examples
#' set.seed(1)
#' pp(cumsum(rnorm(150)), trend = "c")
#' @export
pp <- function(y, trend = "c", bandwidth = "short", level = 0.05, var = NULL) {
  call <- match.call()
  check_level(level)
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  bw <- bandwidth_args(bandwidth)
  bw$rule <- sub("n/100", "(T-1)/100", bw$rule, fixed = TRUE)
  u <- urca::ur.pp(y, type = "Z-tau",
                   model = if (trend == "c") "constant" else "trend",
                   lags = bw$lags, use.lag = bw$use.lag)
  n <- length(y)
  nobs <- n - 1L
  new_ur_test(
    method = "Phillips-Perron test", test = "PP", statistic_name = "Z_tau",
    statistic = unname(u@teststat), null = "unit root",
    alternative = adf_alternatives[[trend]], trend = trend,
    lags = u@lag, bandwidth_rule = bw$rule,
    nobs = nobs, sample = c(2L, n),
    critical_values = tau_critical_values(trend, nobs), tail = "lower",
    calibration = list(
      method = "response_surface", source = mackinnon_source[[trend]],
      note = sprintf("evaluated at N = %d; Z_tau has the Dickey-Fuller tau limit distribution", nobs)),
    level = level, engine = engine_urca(),
    warnings = if (nobs < 25L) "N is below 25." else character(),
    call = call)
}

#' KPSS test of stationarity
#'
#' Kwiatkowski-Phillips-Schmidt-Shin test of the null hypothesis that the
#' series is stationary around a constant level (`trend = "c"`) or a linear
#' trend (`trend = "ct"`), against a unit root. Large values reject.
#' Computed by [urca::ur.kpss()], with a Bartlett kernel.
#'
#' Limitations of this version (backend: urca): bandwidth `"short"`,
#' `"long"` (Schwert's rules) or fixed; no automatic bandwidth; asymptotic
#' critical values.
#'
#' @inheritParams pp
#' @param level Significance level of the reported decision: 0.01, 0.025,
#'   0.05 or 0.10.
#' @return An object of class `ur_test`.
#' @references
#' Kwiatkowski, D., Phillips, P. C. B., Schmidt, P. and Shin, Y. (1992).
#' Testing the null hypothesis of stationarity against the alternative of a
#' unit root. *Journal of Econometrics*, 54, 159-178.
#' @examples
#' set.seed(1)
#' kpss(cumsum(rnorm(150)), trend = "c")
#' @export
kpss <- function(y, trend = "c", bandwidth = "short", level = 0.05, var = NULL) {
  call <- match.call()
  check_level(level, c(0.01, 0.025, 0.05, 0.10))
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  bw <- bandwidth_args(bandwidth)
  bw$rule <- sub("n/100", "T/100", bw$rule, fixed = TRUE)
  u <- urca::ur.kpss(y, type = if (trend == "c") "mu" else "tau",
                     lags = bw$lags, use.lag = bw$use.lag)
  cv <- u@cval[1, ]
  names(cv) <- c(`10pct` = "0.1", `5pct` = "0.05", `2.5pct` = "0.025",
                 `1pct` = "0.01")[colnames(u@cval)]
  n <- length(y)
  stationary <- c(c = "stationary around a constant level",
                  ct = "stationary around a linear trend")[[trend]]
  new_ur_test(
    method = "KPSS test of stationarity", test = "KPSS",
    statistic_name = if (trend == "c") "eta_mu" else "eta_tau",
    statistic = unname(u@teststat), null = stationary,
    alternative = "unit root", trend = trend,
    lags = u@lag, bandwidth_rule = bw$rule,
    nobs = n, sample = c(1L, n),
    critical_values = cv, tail = "upper",
    calibration = list(method = "table",
                       source = "Kwiatkowski et al. (1992), Table 1",
                       note = "asymptotic critical values"),
    level = level, engine = engine_urca(), call = call)
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
#' Limitations of this version (backend: urca): critical values as provided
#' by urca, rounded to two decimals; no p-value. **These critical values are
#' too liberal in finite samples**: at the 5% level the test rejects a true
#' unit root about 8% of the time with \eqn{T = 100} and a constant (10%
#' with \eqn{T = 50}). Finite-sample critical values are planned.
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
  cv <- u@cval[1, ]
  names(cv) <- c(`1pct` = "0.01", `5pct` = "0.05", `10pct` = "0.1")[colnames(u@cval)]
  n <- length(y)
  nobs <- n - k - 1L
  calibration <- if (trend == "c") {
    list(method = "response_surface",
         source = "MacKinnon (1991) Dickey-Fuller surface without deterministic terms, as computed by urca",
         note = "rounded to two decimals by urca")
  } else {
    list(method = "table",
         source = "Elliott, Rothenberg and Stock (1996), Table 1, as computed by urca",
         note = "rounded to two decimals by urca")
  }
  new_ur_test(
    method = "DF-GLS test", test = "DF-GLS",
    statistic_name = "tau_GLS", statistic = unname(u@teststat),
    null = "unit root", alternative = adf_alternatives[[trend]], trend = trend,
    lags = k, lag_selection = selection, nobs = nobs, sample = c(k + 2L, n),
    critical_values = cv, tail = "lower",
    calibration = calibration, level = level, engine = engine_urca(),
    warnings = c(notes, if (nobs < 25L) "N is below 25." else character()),
    call = call)
}
