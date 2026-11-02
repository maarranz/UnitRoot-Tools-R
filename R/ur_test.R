# The result contract (spec: "The result contract")
#
# Every test in the package returns an object of class `ur_test`, whatever
# engine computed it. Tests build it with new_ur_test(); users do not.

calibration_methods <- c("table", "response_surface", "conventional_t", "bootstrap")
required_levels <- c(0.01, 0.05, 0.10)

#' Unit-root test results
#'
#' Objects of class `ur_test` hold the result of a single unit-root or
#' stationarity test in the common format defined by the package
#' specification. They are created by the test functions, not directly
#' by users; `new_ur_test()` is exported so that the format is documented
#' and can be used by extensions.
#'
#' @param method Full name of the test, e.g. `"Augmented Dickey-Fuller test"`.
#' @param test Test family, e.g. `"ADF"`, `"KPSS"`.
#' @param statistic_name Name of the statistic, e.g. `"tau3"`, `"phi3"`.
#' @param statistic Value of the statistic.
#' @param null,alternative Hypotheses in words, e.g. `"unit root"` and
#'   `"stationary around a linear trend"`.
#' @param trend Deterministic specification (canonical code or alias; see
#'   [match_trend()]).
#' @param lags Number of lagged differences (ADF-type tests) or bandwidth
#'   (PP, KPSS); `NA` when not applicable.
#' @param bandwidth_rule For long-run variance based tests, the rule that
#'   set the bandwidth (e.g. `"Schwert short"`); `NULL` otherwise.
#' @param lag_selection `NULL` when lags were fixed by the user, otherwise
#'   a list with elements `method`, `kmax` and, optionally, `table` (the
#'   full criterion table) and `detrending`.
#' @param nobs Effective sample size \eqn{N}.
#' @param sample Integer vector of length 2: first and last \eqn{t} of the
#'   effective sample.
#' @param critical_values Named numeric vector; names are significance
#'   levels (e.g. `"0.01"`, `"0.05"`, `"0.1"`). Levels 0.01, 0.05 and 0.10
#'   are required.
#' @param tail `"lower"` if the test rejects for small values (unit-root
#'   tests), `"upper"` if it rejects for large values (KPSS, \eqn{\phi}).
#' @param p_value P-value, or `NA` when the calibration does not provide
#'   one. Never interpolated from tabulated critical values.
#' @param calibration List with `method` (one of `"table"`,
#'   `"response_surface"`, `"conventional_t"`, `"bootstrap"`), `source`
#'   and, optionally, `note` (e.g. an approximation used).
#' @param level Significance level for the reported decision.
#' @param other_statistics Named numeric vector of accompanying statistics
#'   (e.g. the \eqn{\phi} statistics reported with \eqn{\tau_3}).
#' @param other_critical_values Optional named list, one element per entry
#'   of `other_statistics` that has critical values: a named numeric vector
#'   in the same format as `critical_values`, optionally with a `source`
#'   attribute. These statistics are upper-tail (reject for large values).
#' @param engine Engine that computed the result, e.g. `"internal 0.1.0"`
#'   or `"urca 1.3-4"`.
#' @param warnings Character vector of warnings (small samples,
#'   approximations).
#' @param regression Optional fitted regression (list); may be `NULL`.
#' @param call The call that produced the result.
#'
#' @return An object of class `ur_test`.
#' @export
new_ur_test <- function(method, test, statistic_name, statistic,
                        null, alternative, trend,
                        lags = NA_integer_, bandwidth_rule = NULL,
                        lag_selection = NULL,
                        nobs, sample,
                        critical_values, tail = c("lower", "upper"),
                        p_value = NA_real_,
                        calibration, level = 0.05,
                        other_statistics = numeric(),
                        other_critical_values = list(),
                        engine, warnings = character(),
                        regression = NULL, call = NULL) {
  tail <- match.arg(tail)
  x <- list(
    method = method,
    test = test,
    statistic_name = statistic_name,
    statistic = statistic,
    null = null,
    alternative = alternative,
    trend = match_trend(trend),
    lags = as.integer(lags),
    bandwidth_rule = bandwidth_rule,
    lag_selection = lag_selection,
    nobs = as.integer(nobs),
    sample = as.integer(sample),
    critical_values = critical_values,
    tail = tail,
    p_value = as.numeric(p_value),
    calibration = calibration,
    level = level,
    reject = NA,
    other_statistics = other_statistics,
    other_critical_values = other_critical_values,
    engine = engine,
    warnings = as.character(warnings),
    regression = regression,
    call = call
  )
  class(x) <- "ur_test"
  validate_ur_test(x)
  x$reject <- ur_decision(x$statistic, x$critical_values, x$tail,
                          x$p_value, x$level)
  x
}

#' @rdname new_ur_test
#' @param x An object of class `ur_test`.
#' @return `validate_ur_test()` returns `x` invisibly, or stops with an
#'   informative error.
#' @export
validate_ur_test <- function(x) {
  fail <- function(...) stop("Invalid ur_test: ", ..., call. = FALSE)
  is_string <- function(v) is.character(v) && length(v) == 1L && !is.na(v) && nzchar(v)
  is_number <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)

  for (f in c("method", "test", "statistic_name", "null", "alternative", "engine")) {
    if (!is_string(x[[f]])) fail("`", f, "` must be a non-empty string.")
  }
  if (!is_number(x$statistic)) fail("`statistic` must be a finite number.")

  if (length(x$lags) != 1L || (!is.na(x$lags) && x$lags < 0L)) {
    fail("`lags` must be a non-negative integer or NA.")
  }
  if (!is.null(x$bandwidth_rule) && !is_string(x$bandwidth_rule)) {
    fail("`bandwidth_rule` must be NULL or a string.")
  }
  if (!is.null(x$lag_selection)) {
    ls <- x$lag_selection
    if (!is.list(ls) || !is_string(ls$method) || !is_number(ls$kmax)) {
      fail("`lag_selection` must be NULL or a list with `method` and `kmax`.")
    }
    if (is.na(x$lags) || x$lags > ls$kmax) {
      fail("selected `lags` must lie between 0 and `lag_selection$kmax`.")
    }
  }

  if (length(x$nobs) != 1L || is.na(x$nobs) || x$nobs < 1L) {
    fail("`nobs` must be a positive integer.")
  }
  if (length(x$sample) != 2L || anyNA(x$sample) || x$sample[1] > x$sample[2]) {
    fail("`sample` must be two integers (first, last) with first <= last.")
  }
  if (x$sample[2] - x$sample[1] + 1L != x$nobs) {
    fail("`nobs` must equal the length of `sample`.")
  }

  cv <- x$critical_values
  if (!is.numeric(cv) || is.null(names(cv)) || anyNA(cv) || any(!is.finite(cv))) {
    fail("`critical_values` must be a named numeric vector of finite values.")
  }
  lv <- suppressWarnings(as.numeric(names(cv)))
  if (anyNA(lv) || any(lv <= 0 | lv >= 1)) {
    fail("names of `critical_values` must be significance levels in (0, 1).")
  }
  if (!all(vapply(required_levels, function(l) any(abs(lv - l) < 1e-12), logical(1)))) {
    fail("`critical_values` must include the levels 0.01, 0.05 and 0.10.")
  }
  ord <- order(lv)
  d <- diff(cv[ord])
  ok <- if (x$tail == "lower") all(d > 0) else all(d < 0)
  if (!ok) {
    fail("`critical_values` are not ordered consistently with a ",
         x$tail, "-tail test.")
  }

  if (length(x$p_value) != 1L ||
      (!is.na(x$p_value) && (x$p_value < 0 || x$p_value > 1))) {
    fail("`p_value` must be NA or a number in [0, 1].")
  }

  cal <- x$calibration
  if (!is.list(cal) || !is_string(cal$method) ||
      !cal$method %in% calibration_methods || !is_string(cal$source)) {
    fail("`calibration` must be a list with `method` (one of ",
         paste0('"', calibration_methods, '"', collapse = ", "),
         ") and `source`.")
  }
  if (!is.null(cal$note) && !is_string(cal$note)) {
    fail("`calibration$note` must be NULL or a string.")
  }

  if (!is_number(x$level) || x$level <= 0 || x$level >= 1) {
    fail("`level` must be a number in (0, 1).")
  }
  if (!any(abs(lv - x$level) < 1e-12) && is.na(x$p_value)) {
    fail("no critical value at `level` = ", x$level,
         " and no p-value: the decision cannot be made.")
  }

  os <- x$other_statistics
  if (!is.numeric(os) || (length(os) > 0L && is.null(names(os)))) {
    fail("`other_statistics` must be a named numeric vector.")
  }
  ocv <- x$other_critical_values
  if (!is.list(ocv) || (length(ocv) > 0L && (is.null(names(ocv)) ||
      !all(names(ocv) %in% names(os))))) {
    fail("`other_critical_values` must be a list named after entries of `other_statistics`.")
  }
  for (nm in names(ocv)) {
    v <- ocv[[nm]]
    lvv <- suppressWarnings(as.numeric(names(v)))
    if (!is.numeric(v) || anyNA(lvv) || any(!is.finite(v))) {
      fail("`other_critical_values$", nm, "` must be a named numeric vector of finite values.")
    }
  }
  if (!is.null(x$regression) && !is.list(x$regression)) {
    fail("`regression` must be NULL or a list.")
  }
  invisible(x)
}

# Decision at `level`: from the critical value at that level when tabulated,
# otherwise from the p-value.
ur_decision <- function(statistic, critical_values, tail, p_value, level) {
  lv <- as.numeric(names(critical_values))
  i <- which(abs(lv - level) < 1e-12)
  if (length(i) == 1L) {
    cv <- critical_values[[i]]
    if (tail == "lower") statistic < cv else statistic > cv
  } else {
    p_value < level
  }
}

# Engine tag for results computed by this package.
engine_internal <- function() {
  paste("internal", unname(getNamespaceVersion("unitroottools")))
}

# ---- printing ---------------------------------------------------------------

fmt_num <- function(v, digits = 3) formatC(v, format = "f", digits = digits)

fmt_level <- function(l) paste0(format(100 * l, trim = TRUE), "%")

fmt_pvalue <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) "< 0.001" else fmt_num(p, 3)
}

fmt_cvs <- function(cv) {
  lv <- as.numeric(names(cv))
  o <- order(lv)
  paste(sprintf("%s: %s", vapply(lv[o], fmt_level, ""), fmt_num(cv[o])),
        collapse = "   ")
}

lags_line <- function(x) {
  if (is.na(x$lags)) return(NULL)
  if (!is.null(x$bandwidth_rule)) {
    return(sprintf("Bandwidth: %d (%s)", x$lags, x$bandwidth_rule))
  }
  if (!is.null(x$lag_selection)) {
    return(sprintf("Lags: %d (%s, kmax = %d)", x$lags,
                   x$lag_selection$method, as.integer(x$lag_selection$kmax)))
  }
  sprintf("Lags: %d (fixed)", x$lags)
}

decision_line <- function(x) {
  sprintf("Decision at %s: %s H0 (%s)", fmt_level(x$level),
          if (isTRUE(x$reject)) "reject" else "do not reject", x$null)
}

#' @export
print.ur_test <- function(x, ...) {
  cat(sprintf("%s (%s)\n\n", x$method, x$statistic_name))
  cat(sprintf("  H0: %s\n  H1: %s\n", x$null, x$alternative))
  cat(sprintf("  Deterministics: %s\n", trend_labels[[x$trend]]))
  ll <- lags_line(x)
  if (!is.null(ll)) cat("  ", ll, "\n", sep = "")
  cat(sprintf("  Effective sample: N = %d (t = %d, ..., %d)\n\n",
              x$nobs, x$sample[1], x$sample[2]))
  cat(sprintf("  Statistic:        %s\n", fmt_num(x$statistic)))
  cat(sprintf("  Critical values:  %s\n", fmt_cvs(x$critical_values)))
  if (!is.na(x$p_value)) cat(sprintf("  p-value:          %s\n", fmt_pvalue(x$p_value)))
  cat("\n  ", decision_line(x), "\n", sep = "")
  for (w in x$warnings) cat("  Warning: ", w, "\n", sep = "")
  invisible(x)
}

#' @export
summary.ur_test <- function(object, ...) {
  structure(list(result = object), class = "summary.ur_test")
}

#' @export
print.summary.ur_test <- function(x, ...) {
  r <- x$result
  print(r)
  cal <- r$calibration
  cat("\nCalibration\n")
  cat(sprintf("  Method: %s\n  Source: %s\n", cal$method, cal$source))
  if (!is.null(cal$note)) cat(sprintf("  Note:   %s\n", cal$note))
  cat(sprintf("  p-value: %s\n",
              if (is.na(r$p_value)) "not provided by this calibration"
              else fmt_pvalue(r$p_value)))
  if (length(r$other_statistics) > 0L) {
    cat("\nAccompanying statistics\n")
    for (n in names(r$other_statistics)) {
      cv <- r$other_critical_values[[n]]
      if (is.null(cv)) {
        cat(sprintf("  %-6s %9s\n", n, fmt_num(r$other_statistics[[n]])))
      } else {
        lv <- as.numeric(names(cv))
        i5 <- which(abs(lv - r$level) < 1e-12)
        dec <- if (length(i5) == 1L) {
          if (r$other_statistics[[n]] > cv[[i5]]) "  reject" else "  do not reject"
        } else ""
        cat(sprintf("  %-6s %9s   %s%s\n", n, fmt_num(r$other_statistics[[n]]),
                    fmt_cvs(cv), dec))
        if (!is.null(attr(cv, "source"))) {
          cat(sprintf("  %-6s %9s   [%s]\n", "", "", attr(cv, "source")))
        }
      }
    }
  }
  ls <- r$lag_selection
  if (!is.null(ls)) {
    cat(sprintf("\nLag selection: %s, kmax = %d", ls$method, as.integer(ls$kmax)))
    if (!is.null(ls$detrending)) cat(sprintf(", detrending: %s", ls$detrending))
    if (!is.null(ls$final_sample)) cat(sprintf(", test regression on the %s sample", ls$final_sample))
    cat("\n")
    if (!is.null(ls$table)) {
      tab <- ls$table
      shown <- tab
      for (cn in setdiff(names(tab), "k")) shown[[cn]] <- sprintf("%.4f", tab[[cn]])
      shown[[" "]] <- ifelse(tab$k == r$lags, "  <-", "")
      print(shown, row.names = FALSE, right = TRUE)
    }
  }
  if (!is.null(r$regression$coefficients)) {
    cat("\nTest regression\n")
    print(r$regression$coefficients, digits = 4)
  }
  cat(sprintf("\nEngine: %s\n", r$engine))
  invisible(x)
}
