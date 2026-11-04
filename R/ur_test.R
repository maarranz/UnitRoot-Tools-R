# The result contract (spec: "The result contract")
#
# Every test in the package returns an object of class `ur_test`, whatever
# engine computed it. Tests build it with new_ur_test(); users do not.

calibration_methods <- c("table", "response_surface", "limit_distribution",
                         "conventional_t", "bootstrap")
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
#'   `"response_surface"`, `"limit_distribution"`, `"conventional_t"`,
#'   `"bootstrap"`), `source`
#'   and, optionally, `note` (e.g. an approximation used) and
#'   `p_value_source` (when the p-value comes from a different source than
#'   the critical values).
#' @param level Significance level for the reported decision.
#' @param other_statistics Named numeric vector of accompanying statistics
#'   (e.g. the \eqn{\phi} statistics reported with \eqn{\tau_3}).
#' @param other_critical_values Optional named list, one element per entry
#'   of `other_statistics` that has critical values: a named numeric vector
#'   in the same format as `critical_values`, optionally with attributes
#'   `source`, `tail` (`"upper"`, the default, or `"lower"`) and `p_value`.
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
    lags = as_whole(lags, "lags", allow_na = TRUE),
    bandwidth_rule = bandwidth_rule,
    lag_selection = lag_selection,
    nobs = as_whole(nobs, "nobs"),
    sample = as_whole(sample, "sample", length = 2L),
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
  validate_fields(x)
  x$reject <- ur_decision(x$statistic, x$critical_values, x$tail,
                          x$p_value, x$level)
  validate_ur_test(x)
  x   # returned visibly, so that adf(y) etc. print at the console
}

# Whole numbers are converted to integer only if they are exactly whole;
# 2.5 is an error, never 2.
as_whole <- function(v, name, allow_na = FALSE, length = 1L) {
  ok <- is.numeric(v) && !is.logical(v) && length(v) == length &&
    (if (allow_na) all(is.na(v) | (is.finite(v) & v == round(v) & abs(v) <= .Machine$integer.max))
     else all(!is.na(v) & is.finite(v) & v == round(v) & abs(v) <= .Machine$integer.max))
  if (!ok) {
    stop("Invalid ur_test: `", name, "` must be ",
         if (length == 1L) "a whole number" else paste(length, "whole numbers"),
         if (allow_na) " or NA" else "", " (within the integer range).", call. = FALSE)
  }
  as.integer(v)
}

#' @rdname new_ur_test
#' @param x An object of class `ur_test`.
#' @return `validate_ur_test()` returns `x` invisibly, or stops with an
#'   informative error. It checks every field of the contract, including that
#'   the stored decision agrees with the statistic, critical values, tail,
#'   p-value and level.
#' @export
validate_ur_test <- function(x) {
  validate_fields(x)
  if (!is.logical(x$reject) || length(x$reject) != 1L || is.na(x$reject)) {
    stop("Invalid ur_test: `reject` must be TRUE or FALSE.", call. = FALSE)
  }
  expected <- ur_decision(x$statistic, x$critical_values, x$tail,
                          x$p_value, x$level)
  if (!identical(x$reject, expected)) {
    stop("Invalid ur_test: the stored decision (reject = ", x$reject,
         ") does not agree with the statistic, critical values and level.",
         call. = FALSE)
  }
  invisible(x)
}

# Every field except the decision (the constructor computes the decision
# from these fields, then runs the full validator).
validate_fields <- function(x) {
  fail <- function(...) stop("Invalid ur_test: ", ..., call. = FALSE)
  is_string <- function(v) is.character(v) && length(v) == 1L && !is.na(v) && nzchar(v)
  is_number <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  is_int <- function(v, n = 1L) is.integer(v) && length(v) == n

  if (!inherits(x, "ur_test") || !is.list(x)) fail("not an object of class `ur_test`.")
  for (f in c("method", "test", "statistic_name", "null", "alternative", "engine")) {
    if (!is_string(x[[f]])) fail("`", f, "` must be a non-empty string.")
  }
  if (!is_number(x$statistic)) fail("`statistic` must be a finite number.")
  if (!is_string(x$trend) || !x$trend %in% c("n", "c", "ct")) {
    fail("`trend` must be one of \"n\", \"c\", \"ct\" (canonical codes).")
  }
  if (!is_string(x$tail) || !x$tail %in% c("lower", "upper")) {
    fail("`tail` must be \"lower\" or \"upper\".")
  }

  if (!is_int(x$lags) || (!is.na(x$lags) && x$lags < 0L)) {
    fail("`lags` must be a non-negative integer or NA.")
  }
  if (!is.null(x$bandwidth_rule) && !is_string(x$bandwidth_rule)) {
    fail("`bandwidth_rule` must be NULL or a string.")
  }
  if (!is.null(x$lag_selection)) {
    ls <- x$lag_selection
    if (!is.list(ls) || !is_string(ls$method) || !is_number(ls$kmax) ||
        ls$kmax < 0 || ls$kmax != round(ls$kmax)) {
      fail("`lag_selection` must be NULL or a list with `method` and a whole-number `kmax`.")
    }
    if (is.na(x$lags) || x$lags > ls$kmax) {
      fail("selected `lags` must lie between 0 and `lag_selection$kmax`.")
    }
  }

  if (!is_int(x$nobs) || is.na(x$nobs) || x$nobs < 1L) {
    fail("`nobs` must be a positive integer.")
  }
  if (!is_int(x$sample, 2L) || anyNA(x$sample) || x$sample[1] < 1L ||
      x$sample[1] > x$sample[2]) {
    fail("`sample` must be two integers (first, last) with 1 <= first <= last.")
  }
  if (x$sample[2] - x$sample[1] + 1L != x$nobs) {
    fail("`nobs` must equal the length of `sample`.")
  }

  check_cv <- function(cv, what, tail) {
    if (!is.numeric(cv) || is.null(names(cv)) || anyNA(cv) || any(!is.finite(cv))) {
      fail("`", what, "` must be a named numeric vector of finite values.")
    }
    lv <- suppressWarnings(as.numeric(names(cv)))
    if (anyNA(lv) || any(lv <= 0 | lv >= 1) || anyDuplicated(lv)) {
      fail("names of `", what, "` must be distinct significance levels in (0, 1).")
    }
    if (!all(vapply(required_levels, function(l) any(abs(lv - l) < 1e-12), logical(1)))) {
      fail("`", what, "` must include the levels 0.01, 0.05 and 0.10.")
    }
    d <- diff(cv[order(lv)])
    ok <- if (tail == "lower") all(d > 0) else all(d < 0)
    if (!ok) {
      fail("`", what, "` are not ordered consistently with a ", tail, "-tail test.")
    }
    lv
  }
  # A bootstrap calibration may give only a p-value (no critical values).
  bootstrap_only <- is.list(x$calibration) && identical(x$calibration$method, "bootstrap") &&
    is.numeric(x$critical_values) && length(x$critical_values) == 0L
  lv <- if (bootstrap_only) {
    if (!is_number(x$p_value)) {
      fail("a bootstrap result without critical values must have a p-value.")
    }
    numeric()
  } else {
    check_cv(x$critical_values, "critical_values", x$tail)
  }

  if (!is.numeric(x$p_value) || length(x$p_value) != 1L ||
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
  if (!is.null(cal$p_value_source) && !is_string(cal$p_value_source)) {
    fail("`calibration$p_value_source` must be NULL or a string.")
  }

  if (!is_number(x$level) || x$level <= 0 || x$level >= 1) {
    fail("`level` must be a number in (0, 1).")
  }
  if (!any(abs(lv - x$level) < 1e-12) && is.na(x$p_value)) {
    fail("no critical value at `level` = ", x$level,
         " and no p-value: the decision cannot be made.")
  }

  os <- x$other_statistics
  if (!is.numeric(os) || any(!is.finite(os)) ||
      (length(os) > 0L && (is.null(names(os)) || any(!nzchar(names(os))) ||
                           anyDuplicated(names(os))))) {
    fail("`other_statistics` must be a named numeric vector of finite values.")
  }
  ocv <- x$other_critical_values
  if (!is.list(ocv) || (length(ocv) > 0L && (is.null(names(ocv)) ||
      !all(names(ocv) %in% names(os))))) {
    fail("`other_critical_values` must be a list named after entries of `other_statistics`.")
  }
  for (nm in names(ocv)) {
    tl <- other_tail(ocv[[nm]])
    if (!tl %in% c("lower", "upper")) {
      fail("`other_critical_values$", nm, "` has an invalid `tail` attribute.")
    }
    check_cv(ocv[[nm]], paste0("other_critical_values$", nm), tl)
    pv <- attr(ocv[[nm]], "p_value")
    if (!is.null(pv) && (!is_number(pv) || pv < 0 || pv > 1)) {
      fail("the `p_value` attribute of `other_critical_values$", nm,
           "` must be a number in [0, 1].")
    }
  }

  if (!is.character(x$warnings) || anyNA(x$warnings)) {
    fail("`warnings` must be a character vector.")
  }
  if (!is.null(x$regression) && !is.list(x$regression)) {
    fail("`regression` must be NULL or a list.")
  }
  if (!is.null(x$call) && !is.call(x$call)) {
    fail("`call` must be NULL or a call.")
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
  if (p <= 0.001) "< 0.001" else if (p >= 0.999) "> 0.999" else fmt_num(p, 3)
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
    return(sprintf("Long-run variance: %s", x$bandwidth_rule))
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
  cat(sprintf("  Critical values:  %s\n",
              if (length(x$critical_values)) fmt_cvs(x$critical_values)
              else "none (the decision uses the bootstrap p-value)"))
  cat(sprintf("  Source:           %s\n", x$calibration$source))
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
              else if (is.null(cal$p_value_source)) fmt_pvalue(r$p_value)
              else sprintf("%s [%s]", fmt_pvalue(r$p_value), cal$p_value_source)))
  if (length(r$other_statistics) > 0L) {
    cat("\nAccompanying statistics\n")
    for (n in names(r$other_statistics)) {
      cv <- r$other_critical_values[[n]]
      if (is.null(cv)) {
        cat(sprintf("  %-7s %9s\n", n, fmt_num(r$other_statistics[[n]])))
      } else {
        lv <- as.numeric(names(cv))
        i5 <- which(abs(lv - r$level) < 1e-12)
        dec <- if (length(i5) == 1L) {
          if (other_reject(r$other_statistics[[n]], cv, i5)) "  reject" else "  do not reject"
        } else ""
        pv <- attr(cv, "p_value")
        cat(sprintf("  %-7s %9s   %s%s%s\n", n, fmt_num(r$other_statistics[[n]]),
                    fmt_cvs(cv), if (is.null(pv)) "" else paste0("   p = ", fmt_pvalue(pv)), dec))
        if (!is.null(attr(cv, "source"))) {
          cat(sprintf("  %-7s %9s   [%s]\n", "", "", attr(cv, "source")))
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

# Tail of an accompanying statistic (attribute of its critical values).
other_tail <- function(cv) {
  tl <- attr(cv, "tail")
  if (is.null(tl)) "upper" else tl
}

other_reject <- function(stat, cv, i) {
  if (other_tail(cv) == "lower") stat < cv[[i]] else stat > cv[[i]]
}
