# Recursive and rolling ADF statistics (spec: Decision D20)

#' Recursive and rolling ADF tests
#'
#' The ADF statistic computed on a sequence of subsamples, to see how the
#' evidence about a unit root changes over time.
#'
#' * `window = "recursive"`: fixed start, moving end; observations
#'   \eqn{1, \dots, e} for \eqn{e = n_0, \dots, T}.
#' * `window = "rolling"`: windows of fixed width \eqn{w}; observations
#'   \eqn{e - w + 1, \dots, e}.
#'
#' Default sizes depend on the frequency \eqn{s} of the data (observations
#' per year): \eqn{n_0 = \max(5s, 30)} and \eqn{w = \max(10s, 40)}
#' observations, i.e. 60 and 120 for monthly data, 30 and 40 for quarterly
#' and annual data. The floors keep each regression reasonably sized.
#'
#' **Lags.** By default the number of lagged differences is selected once,
#' on the full sample (as in [adf()]), and kept fixed in every window, so that
#' the statistic moves only because the data change. Note that this uses the
#' whole sample: the lag order of an early window depends on later
#' observations (look-ahead). With `reselect = TRUE` it is selected again in
#' each window, using only that window's data (`max_lags`, if given, caps
#' every window).
#'
#' **Critical values are nominal, pointwise references.** Each window's
#' statistic is compared with the MacKinnon critical value for its own sample
#' size, as in [adf()]. This is a reference for a test at a single date chosen
#' in advance, and only approximately even then: the critical values are for
#' a fixed lag order, so they do not account for the full-sample lag choice
#' (look-ahead) or for re-selection in each window. They are not valid for
#' the question "does the statistic cross the line at some date?": the
#' minimum of the sequence has a different distribution, with larger
#' critical values in absolute value (Banerjee, Lumsdaine and Stock 1992), so
#' a few crossings are expected even under a unit root. Read the plot as a
#' description of how the evidence evolves, not as a test of when it
#' changed.
#'
#' @inheritParams adf
#' @param window `"recursive"` (default) or `"rolling"`.
#' @param min_size First (smallest) sample of the recursive sequence, in
#'   observations; default \eqn{\max(5s, 30)}.
#' @param width Width of the rolling windows, in observations; default
#'   \eqn{\max(10s, 40)}.
#' @param reselect Select the lags again in each window (default `FALSE`).
#' @param frequency Observations per year, overriding the series' own
#'   frequency (needed for a plain vector that is not annual).
#' @return A data frame of class `ur_recursive`, one row per window: `start`,
#'   `end` (time points), `nobs`, `lags`, `statistic`, `critical_value` (at
#'   `level`, pointwise), `p_value` and `reject`; with `print()` and `plot()`
#'   methods.
#' @references
#' Banerjee, A., Lumsdaine, R. L. and Stock, J. H. (1992). Recursive and
#' sequential tests of the unit-root and trend-break hypotheses: theory and
#' international evidence. *Journal of Business & Economic Statistics*, 10,
#' 271-287.
#' @examples
#' set.seed(1)
#' y <- ts(cumsum(rnorm(240)), start = 2000, frequency = 12)
#' rr <- recursive_adf(y, trend = "c")
#' rr
#' plot(rr)
#' plot(recursive_adf(y, trend = "c", window = "rolling"))
#' @export
recursive_adf <- function(y, trend = "c", lags = NULL, max_lags = NULL,
                          criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                          window = c("recursive", "rolling"),
                          min_size = NULL, width = NULL, reselect = FALSE,
                          frequency = NULL, level = 0.05, var = NULL) {
  call <- match.call()
  criterion <- match_criterion(criterion)
  window <- match.arg(window)
  check_level(level)
  check_flag(reselect, "reselect")
  series_nm <- series_name(substitute(y))
  y_input <- y
  y <- check_series(y, var = var)
  series <- series_info(y_input, var, length(y), series_nm)
  trend <- match_trend(trend)
  n <- length(y)
  notes <- character()

  # frequency
  s <- series_frequency(frequency, series)
  if (is.na(s)) {
    s <- 1L
    if (is.null(min_size) && is.null(width)) {
      notes <- c(notes, paste("frequency unknown: window sizes for annual data",
                              "(pass `frequency =` or a ts object)"))
    }
  }

  # window sizes
  size <- if (window == "recursive") {
    if (is.null(min_size)) as.integer(max(5 * s, 30)) else check_order(min_size, "min_size")
  } else {
    if (is.null(width)) as.integer(max(10 * s, 40)) else check_order(width, "width")
  }
  if (size > n) {
    stop(sprintf("The %s (%d observations) is longer than the series (T = %d).",
                 if (window == "recursive") "first sample" else "window width", size, n),
         call. = FALSE)
  }

  # lags: fixed, or selected once on the full sample
  if (!is.null(lags)) {
    if (!is.null(max_lags)) {
      stop("Give either `lags` (fixed order) or `max_lags` (selection), not both.",
           call. = FALSE)
    }
    k_full <- check_order(lags, "lags")
    if (reselect) stop("`reselect = TRUE` needs lag selection: drop `lags`.", call. = FALSE)
    lag_choice <- "fixed"
    lag_mode <- "fixed"
  } else {
    sel <- lag_table(y, trend, max_lags, criterion)
    k_full <- sel$lags
    lag_mode <- if (reselect) "window" else "full_sample"
    lag_choice <- if (reselect) sprintf("%s, selected in each window (kmax %s)", criterion,
                                        if (is.null(max_lags)) "from the window size"
                                        else sprintf("= %d or the window's largest feasible order",
                                                     check_order(max_lags, "max_lags")))
                  else sprintf("%s on the full sample (kmax = %d), fixed", criterion, sel$kmax)
  }
  p_reg <- c(n = 0L, c = 1L, ct = 2L)[[trend]] + 1L
  if (!reselect && size - k_full - 1L - (p_reg + k_full) < 2L) {
    stop(sprintf(paste("A window of %d observations is too small for the ADF regression",
                       "with %d lags; use a larger %s or fewer lags."),
                 size, k_full, if (window == "recursive") "`min_size`" else "`width`"),
         call. = FALSE)
  }

  ends <- size:n
  starts <- if (window == "recursive") rep(1L, length(ends)) else ends - size + 1L
  lv <- format_level(level)
  rows <- lapply(seq_along(ends), function(i) {
    yi <- y[starts[i]:ends[i]]
    k <- if (!reselect) k_full else {
      # the user's max_lags caps every window (audit of 0.3.0); without it,
      # the default kmax of the window, as for adf() on that subsample
      kmax_w <- if (is.null(max_lags)) NULL
                else min(check_order(max_lags, "max_lags"), max(0L, kmax_cap(length(yi), trend)))
      lag_table(yi, trend, kmax_w, criterion)$lags
    }
    r <- adf_regression(yi, trend, k, k + 2L)
    stat <- unname(r$tau)
    cv <- tau_critical_values(trend, r$nobs)[[lv]]
    c(nobs = r$nobs, lags = k, statistic = stat, critical_value = cv,
      p_value = tau_p_value(stat, trend, r$nobs))
  })
  m <- do.call(rbind, rows)
  out <- data.frame(start = series$time[starts], end = series$time[ends],
                    nobs = as.integer(m[, "nobs"]), lags = as.integer(m[, "lags"]),
                    statistic = m[, "statistic"], critical_value = m[, "critical_value"],
                    p_value = m[, "p_value"])
  out$reject <- out$statistic < out$critical_value
  structure(out, class = c("ur_recursive", "data.frame"),
            settings = list(trend = trend, window = window, size = size, frequency = s,
                            lag_choice = lag_choice, lag_mode = lag_mode, level = level, series = series,
                            notes = notes, call = call))
}

# Time points for printing: dates as they are; ts times (decimal years) as
# "2004-12" (monthly) or "2004 Q4" (quarterly); annual ts times as years.
fmt_time <- function(v, s = NA) {
  if (inherits(v, c("Date", "POSIXt"))) {
    d <- as.Date(v)
    if (isTRUE(s == 12)) return(format(d, "%Y-%m"))
    if (isTRUE(s == 4)) return(paste0(format(d, "%Y"), " Q", (as.integer(format(d, "%m")) - 1L) %/% 3L + 1L))
    if (isTRUE(s == 1)) return(format(d, "%Y"))
    return(format(v))
  }
  if (isTRUE(s %in% c(4, 12)) && is.numeric(v)) {
    yr <- floor(v + 1e-8)
    per <- round((v - yr) * s) + 1
    return(if (s == 12) sprintf("%d-%02d", as.integer(yr), as.integer(per))
           else sprintf("%d Q%d", as.integer(yr), as.integer(per)))
  }
  if (isTRUE(s == 1) && is.numeric(v) && all(v == round(v))) return(format(as.integer(v)))
  format(v)
}

# "N = 780, 1961-01 to 2025-12 (t = 13, ..., 792)" when the series has time
# points other than 1..T; "N = 780 (t = 13, ..., 792)" otherwise.
sample_text <- function(nobs, sample, series) {
  idx <- sprintf("t = %d, ..., %d", sample[1], sample[2])
  tm <- series$time
  dated <- !is.null(tm) && length(tm) >= sample[2] &&
    !(is.numeric(tm) && !inherits(tm, c("Date", "POSIXt")) && identical(as.numeric(tm), as.numeric(seq_along(tm))))
  if (!dated) return(sprintf("N = %d (%s)", nobs, idx))
  sprintf("N = %d, %s to %s (%s)", nobs, fmt_time(tm[sample[1]], series$frequency),
          fmt_time(tm[sample[2]], series$frequency), idx)
}

recursive_columns <- c("start", "end", "nobs", "lags", "statistic",
                       "critical_value", "p_value", "reject")

# Subsetting keeps the class and the settings when the result still has all
# the columns (e.g. x[x$reject, ]), and returns a plain data frame otherwise
# (e.g. x[, c("end", "statistic")]), so that print() and plot() always work.
#' @export
`[.ur_recursive` <- function(x, ...) {
  st <- attr(x, "settings")
  out <- NextMethod()
  if (is.data.frame(out) && all(recursive_columns %in% names(out)) && !is.null(st)) {
    attr(out, "settings") <- st
    class(out) <- c("ur_recursive", "data.frame")
  } else if (is.data.frame(out)) {
    attr(out, "settings") <- NULL
    class(out) <- "data.frame"
  }
  out
}

#' @export
print.ur_recursive <- function(x, ...) {
  st <- attr(x, "settings")
  # anything that lost its settings or columns (e.g. after other packages'
  # verbs) prints as the data frame it is
  if (is.null(st) || !all(recursive_columns %in% names(x)) || nrow(x) == 0L) {
    print(as.data.frame(unclass_recursive(x)), ...)
    return(invisible(x))
  }
  cat(sprintf("\n%s ADF test (%s)\n", if (st$window == "recursive") "Recursive" else "Rolling",
              trend_labels[[st$trend]]))
  if (!is.null(st$series$name) && !is.na(st$series$name)) cat("  Series:", st$series$name, "\n")
  cat(sprintf("  %s: %d observations; %d windows, ending %s to %s\n",
              if (st$window == "recursive") "First sample" else "Window width",
              st$size, nrow(x), fmt_time(x$end[1], st$series$frequency),
              fmt_time(x$end[nrow(x)], st$series$frequency)))
  cat(sprintf("  Lags: %s\n", st$lag_choice))
  cat(sprintf("  Statistic: min %.3f, max %.3f, last %.3f\n",
              min(x$statistic), max(x$statistic), x$statistic[nrow(x)]))
  nrej <- sum(x$reject)
  cat(sprintf("  Below the pointwise %s critical value in %d of %d windows",
              fmt_level(st$level), nrej, nrow(x)))
  if (nrej) cat(sprintf(" (first ending %s, last ending %s)",
                        fmt_time(x$end[which(x$reject)[1]], st$series$frequency),
                        fmt_time(x$end[utils::tail(which(x$reject), 1)], st$series$frequency)))
  cat("\n")
  for (n in st$notes) cat("  Note:", n, "\n")
  cat("\n  Note: the critical values are nominal pointwise references (one date chosen\n",
      "  in advance, fixed lag order); they are not valid for the minimum over dates\n",
      "  (Banerjee, Lumsdaine and Stock 1992).", if (identical(st$lag_mode, "full_sample"))
        " The lag order uses the full sample (look-ahead).", "\n\n", sep = "")
  invisible(x)
}

unclass_recursive <- function(x) {
  attr(x, "settings") <- NULL
  class(x) <- "data.frame"
  x
}

#' @rdname recursive_adf
#' @param x A result of `recursive_adf()`.
#' @param ... Graphical parameters passed to [plot()].
#' @export
plot.ur_recursive <- function(x, ...) {
  st <- attr(x, "settings")
  if (is.null(st) || !all(recursive_columns %in% names(x)) || nrow(x) == 0L) {
    stop("This object no longer has the settings or columns of a recursive_adf() ",
         "result; plot the columns directly.", call. = FALSE)
  }
  t <- x$end
  rng <- range(c(x$statistic, x$critical_value))
  main <- sprintf("%s ADF statistic (%s), %s",
                  if (st$window == "recursive") "Recursive" else "Rolling",
                  trend_labels[[st$trend]],
                  if (st$window == "recursive") sprintf("first sample %d obs.", st$size)
                  else sprintf("window %d obs.", st$size))
  args <- utils::modifyList(
    list(x = t, y = x$statistic, type = "l", ylim = rng + c(-0.05, 0.22) * diff(rng),
         xlab = "End of the sample", ylab = "ADF statistic", main = main, cex.main = 0.95),
    list(...))
  do.call(graphics::plot, args)
  graphics::lines(t, x$critical_value, lty = 2, col = "firebrick")
  if (any(x$reject)) {
    graphics::points(t[x$reject], x$statistic[x$reject], pch = 16, cex = 0.5, col = "firebrick")
  }
  graphics::legend("top", bty = "n", cex = 0.8, horiz = TRUE,
                   legend = c("ADF statistic",
                              sprintf("pointwise %s critical value", fmt_level(st$level))),
                   lty = c(1, 2), col = c("black", "firebrick"))
  invisible(x)
}
