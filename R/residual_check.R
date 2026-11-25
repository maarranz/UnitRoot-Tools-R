# Residual diagnostics of the test regression (spec: Decision D19)
#
# The ADF-type tests assume that k lagged differences leave white-noise
# residuals. residual_check() reports the correlogram of the residuals (ACF,
# PACF, Ljung-Box Q) and Breusch-Godfrey tests, with a lag grid that depends
# on the frequency of the data (3s + 3 lags: three years plus three).

#' Check the residuals of a unit-root test regression
#'
#' Do the lagged differences of an ADF-type test leave white-noise residuals?
#' `residual_check()` computes, for the residuals \eqn{\hat e_t} of the test
#' regression:
#'
#' * the **correlogram**: autocorrelations (AC), partial autocorrelations
#'   (PAC) and the Ljung-Box statistic \eqn{Q(h)}, as in EViews' correlogram;
#' * the **Breusch-Godfrey** test, F form (the main one) and LM form
#'   (\eqn{N R^2}, EViews' "Obs*R-squared"), at a few orders.
#'
#' Breusch-Godfrey is valid in a regression with lagged dependent variables,
#' such as the ADF regression; Ljung-Box is not exactly justified there, so its
#' p-values (from \eqn{\chi^2(h - k)}, only for \eqn{h > k}) are approximate.
#'
#' **Lags.** With \eqn{s} the frequency of the data (observations per year),
#' the correlogram goes up to \eqn{3s + 3} lags (15 for quarterly data, 39 for
#' monthly, 6 for annual), capped at \eqn{N/4} in short samples; the
#' Breusch-Godfrey orders are \eqn{s} and \eqn{3s + 3} (at most the number of
#' correlogram lags). The frequency is taken from the series (`ts`, zoo, xts,
#' tsibble, or a data frame with a date column). A plain vector has no
#' frequency: it is treated as annual, with a note, unless `frequency` is
#' given.
#'
#' **Which tests.** [adf()], [dfgls()], [ng_perron()] and [ers_pt()] (the
#' regression behind their spectral estimate), the terminal test of
#' [ur_strategy()], and a whole strategy (one check per test regression it
#' used). Not [pp()] or [kpss()]: their regressions leave autocorrelated
#' residuals by design, and the kernel estimate of the long-run variance
#' accounts for it.
#'
#' Increasing the number of lags until the residuals pass is a lag-selection
#' rule of its own, with unknown consequences for the size of the unit-root
#' test, so it is not done automatically: when the residuals fail, try a
#' larger `lags` in the test and compare.
#'
#' @param x A result of [adf()], [dfgls()], [ng_perron()], [ers_pt()] or
#'   [ur_strategy()].
#' @param lags Number of correlogram lags; default \eqn{3s + 3}, capped at
#'   \eqn{N/4}.
#' @param bg_orders Orders of the Breusch-Godfrey test; default \eqn{s} and
#'   \eqn{3s + 3}.
#' @param frequency Observations per year, overriding the series' own
#'   frequency (needed for a plain vector that is not annual).
#' @param ... Not used.
#' @return An object of class `ur_residual_check` (for a strategy, a list of
#'   them, of class `ur_residual_checks`), with `print()` and `plot()`
#'   methods. Its components include `correlogram` (data frame: lag, ac, pac,
#'   q, df, p_value), `breusch_godfrey` (data frame: order, f, df1, df2,
#'   p_value_f, lm, p_value_lm) and `white_noise`: `TRUE` if no
#'   Breusch-Godfrey F test rejects at the test's level, `FALSE` if one does,
#'   and `NA` when no Breusch-Godfrey test could be computed.
#' @references
#' Breusch, T. S. (1978). Testing for autocorrelation in dynamic linear
#' models. *Australian Economic Papers*, 17, 334-355.
#'
#' Godfrey, L. G. (1978). Testing against general autoregressive and moving
#' average error models when the regressors include lagged dependent
#' variables. *Econometrica*, 46, 1293-1301.
#'
#' Ljung, G. M. and Box, G. E. P. (1978). On a measure of lack of fit in time
#' series models. *Biometrika*, 65, 297-303.
#' @examples
#' set.seed(1)
#' y <- ts(cumsum(arima.sim(list(ar = 0.6), n = 240)), start = 2000, frequency = 12)
#' r0 <- adf(y, trend = "c", lags = 0)
#' residual_check(r0)          # autocorrelated: one lag is missing
#' rc <- residual_check(adf(y, trend = "c", lags = 1))
#' rc
#' plot(rc)
#' @export
residual_check <- function(x, lags = NULL, bg_orders = NULL, frequency = NULL, ...) {
  UseMethod("residual_check")
}

#' @rdname residual_check
#' @export
residual_check.ur_test <- function(x, lags = NULL, bg_orders = NULL,
                                   frequency = NULL, ...) {
  if (x$test %in% c("PP", "KPSS")) {
    stop(x$method, ": the test regression leaves autocorrelated residuals by ",
         "design (the kernel estimate of the long-run variance accounts for ",
         "it), so its residuals are not checked.", call. = FALSE)
  }
  g <- x$regression
  if (is.null(g) || is.null(g$residuals) || is.null(g$X)) {
    stop("This result does not contain the residuals of its test regression ",
         "(results of unitroottools 0.3.0 or later do).", call. = FALSE)
  }
  check_residuals(g$residuals, g$X, k = x$lags, level = x$level,
                  series = x$series, lags = lags, bg_orders = bg_orders,
                  frequency = frequency,
                  label = sprintf("%s (%s), k = %d", x$method, trend_labels[[x$trend]], x$lags))
}

# Residual diagnostics of a regression: residuals e (length N), regressor
# matrix X, k lagged differences among the regressors.
check_residuals <- function(e, X, k, level, series, lags, bg_orders, frequency, label) {
  N <- length(e)
  p <- NCOL(X)
  notes <- character()

  # frequency
  s <- series_frequency(frequency, series)
  if (is.na(s)) {
    s <- 1L
    notes <- c(notes, paste("frequency unknown: assumed annual (s = 1); pass",
                            "`frequency =` or a ts object"))
  }

  # correlogram lags
  if (!is.null(lags)) {
    H <- check_order(lags, "lags")
    if (H < 1L || H >= N) stop("`lags` must be between 1 and N - 1 = ", N - 1L, ".", call. = FALSE)
  } else {
    H <- 3L * s + 3L
    cap <- max(1L, N %/% 4L)
    if (H > cap) {
      notes <- c(notes, sprintf("%d lags (3s + 3) capped at N/4 = %d for N = %d", H, cap, N))
      H <- cap
    }
  }
  r <- drop(stats::acf(e, lag.max = H, plot = FALSE, demean = TRUE)$acf)[-1L]
  pr <- drop(stats::pacf(e, lag.max = H, plot = FALSE)$acf)
  h <- seq_len(H)
  Q <- N * (N + 2) * cumsum(r^2 / (N - h))
  df_q <- h - k
  p_q <- ifelse(df_q > 0, stats::pchisq(Q, pmax(df_q, 1), lower.tail = FALSE), NA_real_)
  correlogram <- data.frame(lag = h, ac = r, pac = pr, q = Q,
                            df = ifelse(df_q > 0, df_q, NA_integer_), p_value = p_q)

  # Breusch-Godfrey
  orders <- if (!is.null(bg_orders)) {
    vapply(bg_orders, check_order, integer(1), name = "bg_orders")
  } else {
    unique(pmin(c(s, 3L * s + 3L), H))
  }
  orders <- sort(unique(orders[orders >= 1L]))
  too_big <- orders[N - p - orders < 2L]
  if (length(too_big)) {
    notes <- c(notes, sprintf("Breusch-Godfrey order(s) %s dropped: too few observations",
                              paste(too_big, collapse = ", ")))
    orders <- setdiff(orders, too_big)
  }
  bg <- do.call(rbind, lapply(orders, function(m) breusch_godfrey(e, X, m)))
  if (is.null(bg)) bg <- data.frame(order = integer(), f = numeric(), df1 = integer(),
                                    df2 = integer(), p_value_f = numeric(),
                                    lm = numeric(), p_value_lm = numeric())

  structure(list(label = label, nobs = N, lags_test = k, frequency = s,
                 lags = H, level = level, residuals = e, X = X, series = series,
                 correlogram = correlogram,
                 breusch_godfrey = bg, notes = notes,
                 # TRUE / FALSE from the Breusch-Godfrey tests; NA when none
                 # could be computed (no completed diagnostic: audit of 0.3.0)
                 white_noise = if (nrow(bg)) !any(bg$p_value_f < level) else NA),
            class = "ur_residual_check")
}

pct <- function(level) paste0(format(100 * level), "%")

# Breusch-Godfrey test of order m: regression of e on X and e_{t-1..t-m}, with
# residuals before the start of the sample set to zero (as lmtest::bgtest).
breusch_godfrey <- function(e, X, m) {
  N <- length(e)
  E <- vapply(seq_len(m), function(j) c(rep(0, j), e[seq_len(N - j)]), numeric(N))
  E <- matrix(E, nrow = N)
  Z <- cbind(X, E)
  q <- qr(Z)
  rss_a <- sum(qr.resid(q, e)^2)
  rss_r <- sum(e^2)
  df2 <- N - q$rank
  f <- ((rss_r - rss_a) / m) / (rss_a / df2)
  lm <- N * (1 - rss_a / rss_r)
  data.frame(order = m, f = f, df1 = m, df2 = df2,
             p_value_f = stats::pf(f, m, df2, lower.tail = FALSE),
             lm = lm, p_value_lm = stats::pchisq(lm, m, lower.tail = FALSE))
}

#' @export
print.ur_residual_check <- function(x, max_lags = NULL, ...) {
  cat("\nResidual check:", x$label, "\n")
  cat(sprintf("  N = %d residuals; frequency s = %d; correlogram up to %d lags\n\n",
              x$nobs, x$frequency, x$lags))
  cg <- x$correlogram
  if (!is.null(max_lags)) cg <- cg[cg$lag <= max_lags, , drop = FALSE]
  tab <- data.frame(Lag = cg$lag,
                    AC = sprintf("%6.3f", cg$ac), PAC = sprintf("%6.3f", cg$pac),
                    Q = sprintf("%8.2f", cg$q),
                    `Prob` = ifelse(is.na(cg$p_value), "     -", sprintf("%6.3f", cg$p_value)),
                    check.names = FALSE)
  print(tab, row.names = FALSE, right = TRUE)
  cat(sprintf("  Approximate bands for AC and PAC: +/- %.3f (1.96/sqrt(N)).\n", 1.96 / sqrt(x$nobs)))
  cat(sprintf("  Q p-values from chi-squared(h - %d), approximate in a regression with",
              x$lags_test), "lagged dependent variables.\n\n")
  cat("  Breusch-Godfrey test (residuals before the sample set to zero):\n")
  bg <- x$breusch_godfrey
  if (nrow(bg)) {
    print(data.frame(Order = bg$order,
                     F = sprintf("%7.3f", bg$f),
                     df = sprintf("(%d, %d)", bg$df1, bg$df2),
                     `Prob F` = sprintf("%6.3f", bg$p_value_f),
                     `Obs*R2` = sprintf("%7.3f", bg$lm),
                     `Prob Chi2` = sprintf("%6.3f", bg$p_value_lm),
                     check.names = FALSE), row.names = FALSE, right = TRUE)
  } else {
    cat("  (no order could be computed)\n")
  }
  for (n in x$notes) cat("  Note:", n, "\n")
  cat("\n ", residual_conclusion(x), "\n\n")
  invisible(x)
}

residual_conclusion <- function(x) {
  bg <- x$breusch_godfrey
  if (!nrow(bg)) return(paste("Residual check unavailable: no Breusch-Godfrey test could be",
                             "computed (orders too large for the number of observations)."))
  bad <- bg[bg$p_value_f < x$level, , drop = FALSE]
  if (!nrow(bad)) {
    sprintf("No evidence of residual autocorrelation at the %s level (Breusch-Godfrey F).",
            pct(x$level))
  } else {
    pv <- fmt_pvalue(bad$p_value_f[1])
    sprintf(paste("Residual autocorrelation at the %s level (Breusch-Godfrey F(%d), p %s):",
                  "%d lagged difference%s may not be enough; try a larger `lags`."),
            pct(x$level), bad$order[1], if (grepl("^[<>]", pv)) pv else paste("=", pv),
            x$lags_test, if (x$lags_test == 1L) "" else "s")
  }
}

#' @rdname residual_check
#' @param y Not used (for compatibility with the `plot()` generic).
#' @export
plot.ur_residual_check <- function(x, y, ...) {
  old <- graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old))
  band <- 1.96 / sqrt(x$nobs)
  cg <- x$correlogram
  lim <- max(abs(c(cg$ac, cg$pac)), band) * 1.1
  for (which in c("ac", "pac")) {
    v <- cg[[which]]
    defaults <- list(x = cg$lag, y = v, type = "h", lwd = 3, ylim = c(-lim, lim),
                     xlab = "Lag", ylab = if (which == "ac") "ACF" else "PACF",
                     main = if (which == "ac") paste("Residuals:", x$label) else "",
                     col = ifelse(abs(v) > band, "firebrick", "grey30"), cex.main = 0.9)
    user <- list(...)
    if (which == "pac") user$main <- NULL   # a custom title goes on the top panel only
    do.call(graphics::plot, utils::modifyList(defaults, user))
    graphics::abline(h = 0)
    graphics::abline(h = c(-band, band), lty = 2, col = "steelblue")
  }
  invisible(x)
}

#' @rdname residual_check
#' @export
residual_check.ur_strategy <- function(x, lags = NULL, bg_orders = NULL,
                                       frequency = NULL, ...) {
  checks <- x$residual_checks
  if (is.null(checks)) {
    stop("This strategy result does not contain residual checks ",
         "(results of unitroottools 0.3.0 or later do).", call. = FALSE)
  }
  if (!is.null(lags) || !is.null(bg_orders) || !is.null(frequency)) {
    # recompute with the requested settings from the terminal test's data:
    # the stored checks hold the residuals' summaries only, so use the
    # regressions kept in each check
    checks <- lapply(checks, function(ch) {
      check_residuals(ch$residuals, ch$X, k = ch$lags_test, level = ch$level,
                      series = ch$series, lags = lags, bg_orders = bg_orders,
                      frequency = frequency, label = ch$label)
    })
  }
  structure(checks, class = "ur_residual_checks")
}

#' @export
print.ur_residual_checks <- function(x, ...) {
  cat("\nResidual checks of the test regressions used by the strategy\n")
  for (ch in x) print(ch, ...)
  invisible(x)
}

#' @export
plot.ur_residual_checks <- function(x, y, ...) {
  for (ch in x) plot(ch, ...)
  invisible(x)
}

# Warning for a test result whose regression leaves autocorrelated residuals
# (default settings of residual_check(); D19). Empty when they pass or when
# the check does not apply.
residual_warning <- function(r) {
  g <- r$regression
  if (is.null(g$residuals) || is.null(g$X)) return(character())
  ch <- tryCatch(check_residuals(g$residuals, g$X, k = r$lags, level = r$level, series = r$series,
                                 lags = NULL, bg_orders = NULL, frequency = NULL, label = ""),
                 error = function(e) NULL)
  if (is.null(ch)) return("Residual check not done: see residual_check().")
  if (isTRUE(ch$white_noise)) character()
  else paste(residual_conclusion(ch), "See residual_check().")
}

# Adds the residual warning to a result (public test functions only).
with_residual_warning <- function(r) {
  w <- residual_warning(r)
  if (length(w)) r$warnings <- c(r$warnings, w)
  r
}

# One line for summary(): Breusch-Godfrey F p-values at the default orders.
residual_summary_line <- function(r) {
  g <- r$regression
  if (r$test %in% c("PP", "KPSS") || is.null(g$residuals) || is.null(g$X)) return(NULL)
  ch <- tryCatch(check_residuals(g$residuals, g$X, k = r$lags, level = r$level, series = r$series,
                                 lags = NULL, bg_orders = NULL, frequency = NULL, label = ""),
                 error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  bg <- ch$breusch_godfrey
  if (!nrow(bg)) return(NULL)
  sprintf("Residual autocorrelation (Breusch-Godfrey F): %s%s [details: residual_check()]",
          paste(sprintf("order %d p %s", bg$order,
                        vapply(bg$p_value_f, function(pv) {
                          v <- fmt_pvalue(pv); if (grepl("^[<>]", v)) v else paste("=", v)
                        }, "")),
                collapse = ", "),
          if (any(grepl("frequency unknown", ch$notes))) " (frequency unknown: annual assumed)" else "")
}
