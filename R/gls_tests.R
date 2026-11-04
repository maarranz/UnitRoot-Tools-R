# GLS-based tests (spec: Decision D15): the Ng-Perron M tests and the
# Elliott-Rothenberg-Stock point-optimal test, with lag selection by Perron
# and Qu (2007) and the autoregressive spectral density estimator.

gls_cbar <- c(c = -7, ct = -13.5)

gls_deterministics <- function(n, trend) {
  if (trend == "c") cbind(const = rep(1, n)) else cbind(const = 1, trend = seq_len(n))
}

# Quasi-difference with parameter a: x_1, x_t - a x_{t-1}.
quasi_difference <- function(x, a) {
  x <- as.matrix(x)
  rbind(x[1, , drop = FALSE], x[-1, , drop = FALSE] - a * x[-nrow(x), , drop = FALSE])
}

# GLS detrending (Elliott, Rothenberg and Stock 1996): regress the
# quasi-differenced series on the quasi-differenced deterministic terms with
# a = 1 + cbar/T; returns the detrended series and the residual sum of
# squares S(a), as used by the point-optimal test.
gls_detrend <- function(y, trend, a = 1 + gls_cbar[[trend]] / length(y)) {
  Z <- gls_deterministics(length(y), trend)
  fit <- qr_fit(quasi_difference(Z, a), drop(quasi_difference(y, a)))
  list(yd = drop(y - Z %*% fit$coefficients), ssr = fit$rss, a = a)
}

# Autoregressive spectral density at frequency zero of the GLS-detrended
# series (Ng and Perron 2001, eq. 6; Perron and Ng 1996): from the regression
# of diff(yd)_t on yd_{t-1} and k lagged differences, t = k+2..T,
#   s2_AR = sigma2_k / (1 - sum(beta_j))^2,  sigma2_k = SSR / N.
ar_spectral_density <- function(yd, k) {
  d <- adf_design(yd, "n", k, k + 2L)
  f <- qr_fit(d$X, d$response)
  beta <- if (k > 0) f$coefficients[paste0("dy_lag", seq_len(k))] else numeric()
  sigma2 <- f$rss / length(d$response)
  list(s2_ar = sigma2 / (1 - sum(beta))^2, sigma2 = sigma2, sum_beta = sum(beta),
       nobs = length(d$response))
}

# All GLS-based statistics for a given lag order k.
gls_statistics <- function(y, trend, k) {
  T <- length(y)
  cbar <- gls_cbar[[trend]]
  g <- gls_detrend(y, trend)
  g1 <- gls_detrend(y, trend, a = 1)
  sp <- ar_spectral_density(g$yd, k)
  s2 <- sp$s2_ar
  yd <- g$yd
  kappa <- sum(yd[-T]^2) / T^2
  last <- yd[T]^2 / T
  mza <- (last - s2) / (2 * kappa)
  msb <- sqrt(kappa / s2)
  mpt <- if (trend == "c") (cbar^2 * kappa - cbar * last) / s2
         else (cbar^2 * kappa + (1 - cbar) * last) / s2
  pt <- (g$ssr - g$a * g1$ssr) / s2
  c(MZa = mza, MZt = mza * msb, MSB = msb, MPT = mpt, PT = pt,
    s2_AR = s2)
}

# Lag order for the GLS-based tests: fixed, or selected as recommended by
# Perron and Qu (2007), by select_lags() on the OLS-detrended series (D9).
gls_lag_order <- function(y, trend, lags, max_lags, criterion) {
  if (!is.null(lags)) {
    if (!is.null(max_lags)) {
      stop("Give either `lags` (fixed order) or `max_lags` (selection), not both.",
           call. = FALSE)
    }
    k <- check_order(lags, "lags")
    if (length(y) - k - 1L < k + 3L) {   # two residual degrees of freedom
      stop("`lags` = ", k, " is too large for T = ", length(y),
         " (the test regression needs two residual degrees of freedom).", call. = FALSE)
    }
    return(list(k = k, selection = NULL, notes = character()))
  }
  sel <- lag_table(y, trend, max_lags, criterion)
  list(k = sel$lags, notes = sel$notes,
       selection = list(method = criterion, kmax = sel$kmax,
                        detrending = "OLS (Perron and Qu, 2007); statistics on GLS-detrended series",
                        table = sel$table, final_sample = "maximal", notes = sel$notes))
}

gls_test_result <- function(stat_name, value, others, trend, lo, level, call,
                            method, test, rs, extra_note = NULL) {
  T <- lo$T
  k <- lo$k
  cv <- rs_critical_values(rs, T, k)
  pv <- rs_p_value(rs, unname(value), T, k)
  new_ur_test(
    method = method, test = test,
    statistic_name = stat_name, statistic = unname(value),
    null = "unit root", alternative = adf_alternatives[[trend]], trend = trend,
    lags = k, lag_selection = lo$selection, nobs = T - k - 1L, sample = c(k + 2L, T),
    critical_values = cv, tail = "lower",
    p_value = pv,
    calibration = list(
      method = "response_surface",
      source = "unitroottools GLS simulation tables (own simulation)",
      note = paste0(sprintf("finite-sample quantiles for T = %d and k = %d lagged differences, interpolated in 1/T", T, k),
                    if (is.null(extra_note)) "" else paste0("; ", extra_note),
                    if (is.null(p_value_note(pv))) "" else paste0("; ", p_value_note(pv)))),
    level = level, other_statistics = others,
    engine = engine_internal(),
    warnings = c(lo$notes, rs_design_warnings(rs, T, k), gls_heavy_tail_warning(trend, T, k)),
    call = call)
}

# The GLS tests share the autoregressive spectral density s2_AR, which is
# unstable when the lag order is large for the sample size: the null
# distributions then have very heavy tails (e.g. the 5% quantile of MZa at
# T = 25 with 6 lags is about -800, against -8.6 without lags). The critical
# values are still those of the simulation, but such a test can hardly
# reject. Flagged when the 5% quantile of MZa is more than twice its value
# without lags at the same T.
gls_heavy_tail_warning <- function(trend, T, k) {
  if (k == 0L) return(character())
  rs <- gls_rs$MZa[[trend]]
  j <- which.min(abs(rs$probs - 0.05))
  ratio <- rs_quantiles(rs, T, k)[[j]] / rs_quantiles(rs, T, 0L)[[j]]
  if (ratio <= 2) return(character())
  sprintf(paste("With %d lags at T = %d the null distribution is very heavy-tailed",
                "(the 5%% critical value of MZa is %s times its value without lags):",
                "the critical values are valid, but the test can hardly reject.",
                "Consider fewer lags or a longer sample."),
          k, T, format(signif(ratio, 2)))
}

#' Ng-Perron M tests
#'
#' The modified unit-root tests of Ng and Perron (2001), computed on the
#' series detrended by GLS (\eqn{\bar c = -7} or \eqn{-13.5}, as in
#' [dfgls()]), with the autoregressive estimate of the spectral density at
#' frequency zero,
#' \deqn{s^2_{AR} = \hat\sigma^2_k / (1 - \textstyle\sum_{j=1}^{k} \hat\beta_j)^2,}
#' from the regression of \eqn{\Delta\tilde y_t} on \eqn{\tilde y_{t-1}} and
#' \eqn{k} lagged differences. With \eqn{\kappa = T^{-2}\sum_{t=2}^{T}
#' \tilde y_{t-1}^2}:
#' \deqn{MZ_\alpha = \frac{T^{-1}\tilde y_T^2 - s^2_{AR}}{2\kappa}, \quad
#'   MSB = \sqrt{\kappa / s^2_{AR}}, \quad MZ_t = MZ_\alpha \times MSB,}
#' and \eqn{MP_T = (\bar c^2 \kappa - \bar c\, T^{-1}\tilde y_T^2)/s^2_{AR}}
#' (constant) or \eqn{(\bar c^2 \kappa + (1 - \bar c) T^{-1}\tilde y_T^2)/
#' s^2_{AR}} (constant and trend). All four reject the unit root for small
#' values.
#'
#' The number of lags is chosen as for [dfgls()] (MAIC on the OLS-detrended
#' series, Perron and Qu 2007) unless fixed with `lags`. Critical values and
#' p-values come from the package's own finite-sample simulation over the
#' sample size and the number of lags, interpolated in \eqn{1/T}. In small
#' samples with several lags the critical values can be very far from the
#' asymptotic ones (the spectral estimate is unstable), so the p-value is a
#' better guide than a comparison with published tables. With a trend, the
#' limits differ slightly from Ng and Perron's (2001) Table 1; see the
#' specification (D15).
#'
#' When the lag order is large for the sample size, the spectral estimate is
#' unstable and the null distributions are very heavy-tailed (at \eqn{T = 25}
#' with 6 lags the 5% critical value of \eqn{MZ_\alpha} is about -800). The
#' result then carries a warning: the critical values are valid, but the
#' test can hardly reject.
#'
#' @inheritParams dfgls
#' @param statistic The statistic reported as the test: `"MZt"` (default),
#'   `"MZa"`, `"MSB"` or `"MPT"`. The other three are reported in
#'   `other_statistics`.
#' @return An object of class `ur_test`.
#' @references
#' Ng, S. and Perron, P. (2001). Lag length selection and the construction of
#' unit root tests with good size and power. *Econometrica*, 69, 1519-1554.
#'
#' Perron, P. and Qu, Z. (2007). A simple modification to improve the finite
#' sample properties of Ng and Perron's unit root tests. *Economics Letters*,
#' 94, 12-19.
#' @examples
#' set.seed(1)
#' y <- cumsum(rnorm(150))
#' ng_perron(y, trend = "c")
#' ng_perron(y, trend = "ct", statistic = "MZa")
#' @export
ng_perron <- function(y, trend = "c", lags = NULL, max_lags = NULL,
                      criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                      statistic = c("MZt", "MZa", "MSB", "MPT"),
                      level = 0.05, var = NULL) {
  call <- match.call()
  criterion <- match.arg(criterion)
  statistic <- match.arg(statistic)
  check_level(level)
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  lo <- gls_lag_order(y, trend, lags, max_lags, criterion)
  lo$T <- length(y)
  g <- gls_statistics(y, trend, lo$k)
  m <- c("MZa", "MZt", "MSB", "MPT")
  gls_test_result(statistic, g[[statistic]], g[c(setdiff(m, statistic), "s2_AR")],
                  trend, lo, level, call,
                  method = sprintf("Ng-Perron %s test", statistic), test = "Ng-Perron",
                  rs = gls_rs[[statistic]][[trend]])
}

#' ERS point-optimal test
#'
#' The feasible point-optimal test of Elliott, Rothenberg and Stock (1996),
#' \deqn{P_T = \frac{S(\bar\alpha) - \bar\alpha\, S(1)}{s^2_{AR}},}
#' where \eqn{S(a)} is the residual sum of squares of the GLS regression of
#' the series quasi-differenced with parameter \eqn{a} on the
#' quasi-differenced deterministic terms, \eqn{\bar\alpha = 1 + \bar c/T}
#' (\eqn{\bar c = -7} or \eqn{-13.5}), and \eqn{s^2_{AR}} is the
#' autoregressive spectral density estimate of [ng_perron()]. Small values
#' reject the unit root.
#'
#' Lags as for [dfgls()]; critical values and p-value from the package's own
#' finite-sample simulation.
#'
#' @inheritParams dfgls
#' @return An object of class `ur_test`.
#' @references
#' Elliott, G., Rothenberg, T. J. and Stock, J. H. (1996). Efficient tests for
#' an autoregressive unit root. *Econometrica*, 64, 813-836.
#' @examples
#' set.seed(1)
#' ers_pt(cumsum(rnorm(150)), trend = "c")
#' @export
ers_pt <- function(y, trend = "c", lags = NULL, max_lags = NULL,
                   criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                   level = 0.05, var = NULL) {
  call <- match.call()
  criterion <- match.arg(criterion)
  check_level(level)
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  lo <- gls_lag_order(y, trend, lags, max_lags, criterion)
  lo$T <- length(y)
  g <- gls_statistics(y, trend, lo$k)
  gls_test_result("PT", g[["PT"]], g["s2_AR"], trend, lo, level, call,
                  method = "ERS point-optimal test", test = "ERS-PT",
                  rs = gls_rs[["PT"]][[trend]])
}

#' DF-GLS test
#'
#' Elliott-Rothenberg-Stock DF-GLS test of a unit root: the series is
#' detrended by GLS (local-to-unity, \eqn{\bar c = -7} or \eqn{-13.5}) and a
#' Dickey-Fuller regression without deterministic terms is run on the
#' detrended series.
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
#' The statistic is computed by the package (the same GLS detrending as
#' [ng_perron()] and [ers_pt()]); it agrees with [urca::ur.ers()] to rounding
#' error and with an independent implementation of the published formulas.
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
  lo <- gls_lag_order(y, trend, lags, max_lags, criterion)
  k <- lo$k
  n <- length(y)
  stat <- dfgls_statistic(y, trend, k)
  cv <- dfgls_critical_values(trend, n, k)
  pv <- dfgls_p_value(stat, trend, n, k)
  calibration <- list(
    method = "response_surface",
    source = "unitroottools DF-GLS response surfaces (own simulation)",
    note = paste(c(sprintf("finite-sample quantiles for T = %d and k = %d lagged differences", n, k),
                   p_value_note(pv)), collapse = "; "))
  new_ur_test(
    method = "DF-GLS test", test = "DF-GLS",
    statistic_name = "tau_GLS", statistic = stat,
    null = "unit root", alternative = adf_alternatives[[trend]], trend = trend,
    lags = k, lag_selection = lo$selection, nobs = n - k - 1L, sample = c(k + 2L, n),
    critical_values = cv, tail = "lower",
    p_value = pv,
    calibration = calibration, level = level, engine = engine_internal(),
    warnings = c(lo$notes, dfgls_design_warnings(trend, n, k)),
    call = call)
}

# DF-GLS t statistic: t ratio of yd_{t-1} in the regression of diff(yd)_t on
# yd_{t-1} and k lagged differences, t = k+2..T, without deterministic terms,
# where yd is the GLS-detrended series (Elliott, Rothenberg and Stock 1996).
dfgls_statistic <- function(y, trend, k) {
  yd <- gls_detrend(y, trend)$yd
  d <- adf_design(yd, "n", k, k + 2L)
  unname(qr_fit(d$X, d$response)$t_values[["y_lag1"]])
}
