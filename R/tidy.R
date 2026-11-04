# tidy() and glance() methods (spec: "The result contract", Decision D10)
#
# The generics come from the lightweight 'generics' package and are
# re-exported, so tidy() and glance() work without loading broom (and agree
# with broom when it is loaded). The methods return plain data frames.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy summaries of unit-root tests and strategies
#'
#' `tidy()` returns one row per statistic and `glance()` one row per test or
#' strategy, as plain data frames, ready to be combined across series (e.g.
#' with `rbind()`) or written to a file.
#'
#' For a `ur_test`, `tidy()` has a row for the test statistic and one for
#' each accompanying statistic (e.g. \eqn{\phi_2} and \eqn{\phi_3} with
#' \eqn{\tau_3}), with its critical value at the test's `level`, p-value
#' (`NA` when the calibration provides none), decision and source.
#' `glance()` describes the test: method, deterministic terms, lags and how
#' they were chosen, effective sample, decision, calibration and engine.
#'
#' For a `ur_strategy`, `tidy()` returns the history, one row per step (the
#' same as `as.data.frame()`), and `glance()` the outcome of the procedure.
#'
#' @param x A `ur_test` or `ur_strategy`.
#' @param ... Not used.
#' @return A data frame.
#' @examples
#' set.seed(1)
#' y <- cumsum(0.2 + rnorm(150))
#' r <- adf(y, trend = "ct")
#' tidy(r)
#' glance(r)
#' glance(ur_strategy(y, "original"))
#' @name tidy_methods
NULL

#' @rdname tidy_methods
#' @export
tidy.ur_test <- function(x, ...) {
  main <- data.frame(
    statistic_name = x$statistic_name,
    statistic = x$statistic,
    critical_value = cv_at(x$critical_values, x$level),
    p_value = x$p_value,
    reject = x$reject,
    tail = x$tail,
    source = x$calibration$source,
    stringsAsFactors = FALSE
  )
  os <- x$other_statistics
  if (length(os) == 0L) return(main)
  other <- lapply(names(os), function(nm) {
    cv <- x$other_critical_values[[nm]]
    crit <- if (is.null(cv)) NA_real_ else cv_at(cv, x$level)
    data.frame(
      statistic_name = nm,
      statistic = unname(os[[nm]]),
      critical_value = crit,
      p_value = if (is.null(cv) || is.null(attr(cv, "p_value"))) NA_real_
                else attr(cv, "p_value"),
      reject = if (is.na(crit)) NA
               else if (other_tail(cv) == "lower") unname(os[[nm]]) < crit
               else unname(os[[nm]]) > crit,
      tail = if (is.null(cv)) NA_character_ else other_tail(cv),
      source = if (is.null(cv) || is.null(attr(cv, "source"))) NA_character_
               else attr(cv, "source"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, c(list(main), other))
}

#' @rdname tidy_methods
#' @export
glance.ur_test <- function(x, ...) {
  ls <- x$lag_selection
  data.frame(
    method = x$method,
    test = x$test,
    statistic_name = x$statistic_name,
    statistic = x$statistic,
    p_value = x$p_value,
    level = x$level,
    reject = x$reject,
    trend = x$trend,
    lags = x$lags,
    lag_choice = if (!is.null(x$bandwidth_rule)) x$bandwidth_rule
                 else if (is.null(ls)) "fixed" else ls$method,
    kmax = if (is.null(ls)) NA_integer_ else as.integer(ls$kmax),
    nobs = x$nobs,
    sample_start = x$sample[1],
    sample_end = x$sample[2],
    calibration = x$calibration$source,
    engine = x$engine,
    n_warnings = length(x$warnings),
    stringsAsFactors = FALSE
  )
}

#' @rdname tidy_methods
#' @export
tidy.ur_strategy <- function(x, ...) x$history

#' @rdname tidy_methods
#' @export
glance.ur_strategy <- function(x, ...) {
  ls <- x$lag_selection
  data.frame(
    strategy = x$strategy,
    specification = x$specification,
    unit_root_rejected = x$unit_root_rejected,
    conclusion = x$conclusion,
    terminal_test = x$terminal$statistic_name,
    level = x$settings$level,
    lags = x$lags,
    lag_choice = if (is.null(ls)) "fixed" else ls$method,
    nobs = x$nobs,
    sample_start = x$sample[1],
    sample_end = x$sample[2],
    steps = nrow(x$history),
    n_warnings = length(x$warnings),
    stringsAsFactors = FALSE
  )
}

# Critical value at a level, or NA when that level is not tabulated.
cv_at <- function(cv, level) {
  i <- which(abs(as.numeric(names(cv)) - level) < 1e-12)
  if (length(i) == 1L) unname(cv[[i]]) else NA_real_
}
