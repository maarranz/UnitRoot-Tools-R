# KPSS test of stationarity (spec: Decisions D13 and D14): internal
# implementation with the long-run-variance module, and critical values and
# p-values from the limiting distribution, computed numerically.

# P(eta > q) for eta = sum_k lambda_k Z_k^2, by Imhof (1961). The 200 largest
# eigenvalues are used exactly; the rest contribute their mean (their variance
# is below 1e-7, so the probabilities are unaffected to the precision needed).
kpss_upper_prob <- function(q, lam, K = 200L) {
  if (length(lam) > K) {
    q <- q - sum(lam[-seq_len(K)])
    lam <- lam[seq_len(K)]
  }
  K <- length(lam)
  f <- function(u) {
    m <- length(u)
    th <- 0.5 * .colSums(atan(outer(lam, u)), K, m) - 0.5 * q * u
    rho <- exp(0.25 * .colSums(log1p(outer(lam^2, u^2)), K, m))
    sin(th) / (u * rho)
  }
  p <- 0.5 + stats::integrate(f, 0, Inf, subdivisions = 2000L,
                              rel.tol = 1e-8, abs.tol = 1e-12)$value / pi
  min(max(p, 0), 1)
}

kpss_cv_cache <- new.env(parent = emptyenv())

# Quantiles of the limit where the p-value is bounded (P = 0.999 and 0.001),
# as for the simulated calibrations: beyond them the integral is not
# evaluated (for large statistics it does not converge numerically) and the
# p-value is reported as the bound.
kpss_p_range <- function(trend) {
  key <- paste0("range_", trend)
  if (!is.null(kpss_cv_cache[[key]])) return(kpss_cv_cache[[key]])
  lam <- kpss_limit[[c(c = "mu", ct = "tau")[[trend]]]]
  r <- vapply(c(0.999, 0.001), function(p) {
    stats::uniroot(function(q) kpss_upper_prob(q, lam) - p, c(0.005, 5), tol = 1e-10)$root
  }, numeric(1))
  kpss_cv_cache[[key]] <- r
  r
}

kpss_p_value <- function(stat, trend) {
  r <- kpss_p_range(trend)
  if (stat <= r[1]) return(0.999)
  if (stat >= r[2]) return(0.001)
  kpss_upper_prob(stat, kpss_limit[[c(c = "mu", ct = "tau")[[trend]]]])
}

kpss_critical_values <- function(trend) {
  if (!is.null(kpss_cv_cache[[trend]])) return(kpss_cv_cache[[trend]])
  lam <- kpss_limit[[c(c = "mu", ct = "tau")[[trend]]]]
  cv <- vapply(c(0.10, 0.05, 0.025, 0.01), function(p) {
    stats::uniroot(function(q) kpss_upper_prob(q, lam) - p, c(0.01, 5),
                   tol = 1e-10)$root
  }, numeric(1))
  names(cv) <- c("0.1", "0.05", "0.025", "0.01")
  kpss_cv_cache[[trend]] <- cv
  cv
}

#' KPSS test of stationarity
#'
#' Kwiatkowski-Phillips-Schmidt-Shin test of the null hypothesis that the
#' series is stationary around a constant level (`trend = "c"`) or a linear
#' trend (`trend = "ct"`), against a unit root. Large values reject:
#' \deqn{\eta = \frac{\sum_{t=1}^{T} S_t^2}{T^2 \hat\lambda^2}, \qquad
#'   S_t = \sum_{s \le t} \hat e_s,}
#' where \eqn{\hat e_t} are the residuals of the regression of \eqn{y_t} on
#' the deterministic terms and \eqn{\hat\lambda^2} is their long-run variance.
#'
#' **Long-run variance.** Kernel `"bartlett"` (default, as in Kwiatkowski et
#' al.), `"parzen"` or `"qs"` (Quadratic Spectral). Bandwidth `"short"`
#' (default) or `"long"` (Schwert's rules, \eqn{\ell = \lfloor 4 (T/100)^{1/4}
#' \rfloor} or \eqn{\lfloor 12 (T/100)^{1/4} \rfloor} autocovariances, with
#' Bartlett weights \eqn{1 - j/(\ell + 1)}), a fixed number of
#' autocovariances, or chosen automatically: `"andrews"` (Andrews, 1991, AR(1)
#' plug-in) or `"newey-west"` (Newey and West, 1994). Hobijn, Franses and Ooms
#' (2004) recommend `kernel = "qs", bandwidth = "newey-west"` for this test.
#'
#' **Critical values and p-value** come from the limiting distribution, computed numerically
#' of the statistic, computed numerically (the eigenvalues of the covariance
#' kernel of the Brownian bridge, and Imhof's inversion formula). The values
#' tabulated by Kwiatkowski et al. (1992) were simulated and differ from
#' these in the third decimal (e.g. 0.146 instead of 0.148 at 5% with a
#' trend). As for every KPSS implementation, these are asymptotic: the
#' finite-sample size depends on the bandwidth and on the persistence of
#' the series.
#'
#' @inheritParams pp
#' @param kernel `"bartlett"` (default), `"parzen"` or `"qs"`.
#' @param bandwidth `"short"` (default), `"long"`, `"andrews"`,
#'   `"newey-west"`, or a non-negative integer number of autocovariances.
#' @param level Significance level of the reported decision: 0.01, 0.025,
#'   0.05 or 0.10.
#' @return An object of class `ur_test`.
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. *Econometrica*, 59, 817-858.
#'
#' Hobijn, B., Franses, P. H. and Ooms, M. (2004). Generalizations of the KPSS
#' test for stationarity. *Statistica Neerlandica*, 58, 483-502.
#'
#' Imhof, J. P. (1961). Computing the distribution of quadratic forms in
#' normal variables. *Biometrika*, 48, 419-426.
#'
#' Kwiatkowski, D., Phillips, P. C. B., Schmidt, P. and Shin, Y. (1992).
#' Testing the null hypothesis of stationarity against the alternative of a
#' unit root. *Journal of Econometrics*, 54, 159-178.
#'
#' Newey, W. K. and West, K. D. (1994). Automatic lag selection in covariance
#' matrix estimation. *Review of Economic Studies*, 61, 631-653.
#' @examples
#' set.seed(1)
#' y <- cumsum(rnorm(150))
#' kpss(y, trend = "c")
#' kpss(y, trend = "c", kernel = "qs", bandwidth = "newey-west")
#' @export
kpss <- function(y, trend = "c", kernel = "bartlett", bandwidth = "short",
                 level = 0.05, var = NULL) {
  call <- match.call()
  check_level(level, c(0.01, 0.025, 0.05, 0.10))
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  kernel <- check_kernel(kernel)
  n <- length(y)
  t <- seq_len(n)
  Z <- if (trend == "c") cbind(const = rep(1, n)) else cbind(const = 1, trend = t)
  fit <- qr_fit(Z, y)
  e <- fit$residuals
  if (sqrt(mean(e^2)) <= 1e-10 * max(1, sqrt(mean(y^2)))) {
    stop("The residuals of the regression on the deterministic terms are zero: ",
         "the series is constant", if (trend == "ct") " or an exact linear trend", ". ",
         "KPSS is not defined.", call. = FALSE)
  }
  bw <- resolve_bandwidth(bandwidth, e, kernel, n_rule = n, n_label = "T")
  if (bw$b - 1 >= n) {
    stop("`bandwidth` is too large for T = ", n, ".", call. = FALSE)
  }
  v <- long_run_variance(e, kernel, bw$b)
  if (!(v$lrv > 0)) {
    stop("The estimated long-run variance is not positive; try another kernel or bandwidth.",
         call. = FALSE)
  }
  stat <- sum(cumsum(e)^2) / (n^2 * v$lrv)
  pv <- kpss_p_value(stat, trend)
  stationary <- c(c = "stationary around a constant level",
                  ct = "stationary around a linear trend")[[trend]]
  new_ur_test(
    method = "KPSS test of stationarity", test = "KPSS",
    statistic_name = if (trend == "c") "eta_mu" else "eta_tau",
    statistic = stat, null = stationary,
    alternative = "unit root", trend = trend,
    lags = v$lags_used,
    bandwidth_rule = bandwidth_text(kernel, bw),
    nobs = n, sample = c(1L, n),
    critical_values = kpss_critical_values(trend), tail = "upper",
    p_value = pv,
    calibration = list(
      method = "limit_distribution",
      source = "limiting distribution, computed numerically (kernel eigenvalues, Imhof inversion)",
      note = paste(c("asymptotic; Kwiatkowski et al. (1992) tabulate simulated values that differ in the third decimal",
                     if (pv <= 0.001) "the statistic lies beyond the 0.001 quantile: the p-value is below 0.001 (0.001 is reported)"
                     else if (pv >= 0.999) "the statistic lies below the 0.999 quantile: the p-value is above 0.999 (0.999 is reported)"),
                   collapse = "; ")),
    level = level,
    other_statistics = c(gamma0 = unname(v$gamma0), long_run_variance = unname(v$lrv),
                         bandwidth = unname(v$b)),
    engine = engine_internal(), call = call)
}
