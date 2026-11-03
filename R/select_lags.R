# Lag selection (spec: "Lag selection")

lag_criteria <- c("MAIC", "AIC", "BIC", "MBIC")

default_kmax <- function(n) as.integer(floor(12 * (n / 100)^(1 / 4)))

#' Lag order for the ADF regression
#'
#' Chooses the number of lagged differences in the ADF regression by an
#' information criterion. All candidate orders \eqn{k = 0, \dots, k_{\max}}
#' are estimated on the same (common) sample, \eqn{t = k_{\max} + 2, \dots, T},
#' so that the criteria compare models fitted to the same observations.
#' The full criterion table is returned, not only the selected order.
#'
#' * **AIC** and **BIC** are computed from the ADF regression itself (with
#'   its deterministic terms), with \eqn{\hat\sigma^2_k = RSS_k / N}.
#' * **MAIC** and **MBIC** (Ng and Perron, 2001) are computed from the
#'   regression of \eqn{\Delta\tilde y_t} on \eqn{\tilde y_{t-1}} and
#'   \eqn{k} lagged differences, where \eqn{\tilde y_t} is the series
#'   detrended by OLS under `trend`; the penalty includes
#'   \eqn{\tau_N(k) = \hat\sigma_k^{-2} \hat\beta_0^2 \sum \tilde y_{t-1}^2}.
#'
#' Ties are broken in favour of the smaller order.
#'
#' @param y A univariate time series: numeric vector, `ts`, `zoo`/`xts`, or a
#'   data frame, tibble or tsibble (see `var`).
#' @param trend Deterministic specification: `"n"`, `"c"` or `"ct"`, or an
#'   alias accepted by [match_trend()].
#' @param max_lags Largest order considered, \eqn{k_{\max}}. Default:
#'   \eqn{\lfloor 12 (T/100)^{1/4} \rfloor} (Schwert, 1989).
#' @param criterion One of `"MAIC"` (default), `"AIC"`, `"BIC"`, `"MBIC"`.
#' @param var When `y` has several columns, the name of the column to use.
#'
#' @return An object of class `ur_lags`: a list with the selected order
#'   (`lags`), `criterion`, `kmax`, `trend`, `nobs`, `sample` and `table`, a
#'   data frame with one row per order and one column per criterion.
#' @references
#' Ng, S. and Perron, P. (2001). Lag length selection and the construction of
#' unit root tests with good size and power. *Econometrica*, 69, 1519-1554.
#'
#' Schwert, G. W. (1989). Tests for unit roots: a Monte Carlo investigation.
#' *Journal of Business & Economic Statistics*, 7, 147-159.
#' @examples
#' set.seed(1)
#' y <- cumsum(rnorm(150))
#' select_lags(y, trend = "c")
#' @export
select_lags <- function(y, trend = "c", max_lags = NULL,
                        criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                        var = NULL) {
  criterion <- match.arg(criterion)
  y <- check_series(y, var = var)
  trend <- match_trend(trend)
  lag_table(y, trend, max_lags, criterion)
}

# Internal version on a validated numeric series.
lag_table <- function(y, trend, max_lags, criterion) {
  n <- length(y)
  kmax <- if (is.null(max_lags)) default_kmax(n) else check_order(max_lags, "max_lags")
  first <- kmax + 2L
  nobs <- n - first + 1L
  if (nobs < kmax + 4L) {
    stop("The series is too short for max_lags = ", kmax, " (T = ", n,
         "). Use a smaller `max_lags`.", call. = FALSE)
  }

  ytilde <- ols_detrend(y, trend)
  rows <- lapply(0:kmax, function(k) {
    fit <- adf_regression(y, trend, k, first)$fit
    s2 <- fit$rss / nobs

    d <- adf_design(ytilde, "n", k, first)
    fd <- qr_fit(d$X, d$response)
    s2d <- fd$rss / nobs
    tau_k <- fd$coefficients[["y_lag1"]]^2 * sum(d$X[, "y_lag1"]^2) / s2d

    c(k = k,
      AIC = log(s2) + 2 * k / nobs,
      BIC = log(s2) + k * log(nobs) / nobs,
      MAIC = log(s2d) + 2 * (tau_k + k) / nobs,
      MBIC = log(s2d) + log(nobs) * (tau_k + k) / nobs)
  })
  tab <- as.data.frame(do.call(rbind, rows))
  tab$k <- as.integer(tab$k)
  selected <- tab$k[which.min(tab[[criterion]])]   # which.min: first minimum

  structure(list(lags = selected, criterion = criterion, kmax = kmax,
                 trend = trend, nobs = nobs, sample = c(first, n),
                 detrending = if (trend == "n") "none" else "OLS",
                 table = tab),
            class = "ur_lags")
}

# Residuals from the OLS regression of y on the deterministic terms over the
# full sample t = 1, ..., T.
ols_detrend <- function(y, trend) {
  if (trend == "n") return(y)
  t <- seq_along(y)
  Z <- if (trend == "c") cbind(rep(1, length(y))) else cbind(1, t)
  qr.resid(qr(Z), y)
}

#' @export
print.ur_lags <- function(x, ...) {
  cat(sprintf("Lag selection for the ADF regression (%s)\n\n", x$criterion))
  cat(sprintf("  Deterministics: %s\n", trend_labels[[x$trend]]))
  cat(sprintf("  Common sample:  N = %d (t = %d, ..., %d), kmax = %d\n",
              x$nobs, x$sample[1], x$sample[2], x$kmax))
  cat(sprintf("  MAIC/MBIC detrending: %s\n\n", x$detrending))
  tab <- x$table
  mark <- ifelse(tab$k == x$lags, "  <-", "")
  shown <- data.frame(k = tab$k,
                      AIC = sprintf("%.4f", tab$AIC), BIC = sprintf("%.4f", tab$BIC),
                      MAIC = sprintf("%.4f", tab$MAIC), MBIC = sprintf("%.4f", tab$MBIC),
                      ` ` = mark, check.names = FALSE)
  print(shown, row.names = FALSE, right = TRUE)
  cat(sprintf("\n  Selected by %s: k = %d\n", x$criterion, x$lags))
  invisible(x)
}
