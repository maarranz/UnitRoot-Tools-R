# Union-of-rejections unit-root test (Decision D22): Harvey, Leybourne and
# Taylor (2012), with the package's own finite-sample critical values.

# The four component tests, in the order of the simulation kernel
# (union_null_stats() in data-raw/tau_kernel.cpp).
union_components <- c("ADF_c", "ADF_ct", "GLS_c", "GLS_ct")
union_component_labels <- c(
  ADF_c  = "ADF (OLS), constant",
  ADF_ct = "ADF (OLS), constant and trend",
  GLS_c  = "DF-GLS, constant",
  GLS_ct = "DF-GLS, constant and trend")

# Response surface of each component: tau2, tau3 (D21) and DF-GLS (D11), all
# indexed by the length of the series T and the lag order k (maximal sample).
union_component_rs <- function(name) {
  switch(name, ADF_c = tau_rs$c, ADF_ct = tau_rs$ct,
         GLS_c = dfgls_rs$c, GLS_ct = dfgls_rs$ct)
}

# 5% critical values of the components for (T, k): the normalisation of the
# union statistic (also used by data-raw/union_simulate.R).
union_cv5 <- function(T, k) {
  vapply(union_components, function(nm) {
    unname(rs_critical_values(union_component_rs(nm), T, k, levels = 0.05))
  }, numeric(1))
}

union_source <- "unitroottools union response surface (own simulation)"

#' Union-of-rejections unit-root test
#'
#' Tests for a unit root when it is not known whether a linear trend belongs
#' in the model, nor how large the initial condition is: the null is rejected
#' when any of four tests rejects, ADF (OLS detrending) and DF-GLS (GLS
#' detrending), each with a constant and with a constant and trend (Harvey,
#' Leybourne and Taylor, 2012).
#'
#' **Why a union.** OLS-based tests (ADF) keep their power when the initial
#' condition is large, GLS-based tests (DF-GLS) are more powerful when it is
#' small; tests with a trend are valid whether or not the trend is present
#' but lose power when it is not. Choosing one test after looking at the data
#' (or running a sequence of tests, as [ur_strategy()] does) distorts the size;
#' the union rejects when any test rejects, with a critical value enlarged so
#' that the union as a whole has the nominal size.
#'
#' **Statistic.** Each component statistic is divided by the absolute value
#' of its own 5% critical value for the sample size and lag order, so that
#' each rejects at 5% below \eqn{-1}; the union statistic is the smallest,
#' \deqn{U = \min_i \; t_i / |c_i(0.05; T, k)|,}
#' and the test rejects for small values of \eqn{U}. Its critical values and
#' p-value come from the package's own finite-sample response surface for
#' \eqn{U} under a driftless random walk, in the sample size and the number
#' of lags (Decision D22); the papers use asymptotic scaling constants
#' instead. At 5% the union critical value is about \eqn{-1.15}: the union
#' rejects when a component is 15% beyond its own 5% critical value.
#'
#' **Lags.** One lag order \eqn{k} is used by the four tests (the simulation
#' that calibrates \eqn{U} uses a common \eqn{k}). By default it is selected as
#' in [dfgls()], by `criterion` on the ADF regression with constant and trend
#' (MAIC: OLS-detrended series, Perron and Qu, 2007); every regression then
#' uses the maximal sample \eqn{t = k + 2, \dots, T}.
#'
#' **Reading the result.** `other_statistics` holds the four component
#' statistics, each with its own critical values and p-value (shown by
#' `summary()`); the note says which component attains the minimum. The
#' residual checks ([residual_check()]) use the ADF regression with constant
#' and trend.
#'
#' @inheritParams adf
#' @return An object of class `ur_test` (see [new_ur_test()]) with
#'   `test = "Union"` and `trend = "ct"` (the most general specification
#'   among the components).
#' @references
#' Harvey, D. I., Leybourne, S. J. and Taylor, A. M. R. (2009). Unit root
#' testing in practice: dealing with uncertainty over the trend and initial
#' condition. *Econometric Theory*, 25, 587-636.
#'
#' Harvey, D. I., Leybourne, S. J. and Taylor, A. M. R. (2012). Testing for
#' unit roots in the presence of uncertainty over both the trend and the
#' initial condition. *Journal of Econometrics*, 169, 188-195.
#'
#' Smeekes, S. and Taylor, A. M. R. (2012). Bootstrap union tests for unit
#' roots in the presence of nonstationary volatility. *Econometric Theory*,
#' 28, 422-456.
#' @seealso [boot_ur()] with `union = TRUE` for the bootstrap version
#'   (bootUR); [ur_strategy()] for the sequential alternative.
#' @examples
#' set.seed(1)
#' y <- cumsum(rnorm(200))
#' union_test(y)
#' summary(union_test(y, lags = 2))
#' @export
union_test <- function(y, lags = NULL, max_lags = NULL,
                       criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                       level = 0.05, var = NULL) {
  call <- match.call()
  criterion <- match_criterion(criterion)
  check_level(level)
  series_nm <- series_name(substitute(y))
  y_input <- y
  y <- check_series(y, var = var)
  series <- series_info(y_input, var, length(y), series_nm)
  lo <- gls_lag_order(y, "ct", lags, max_lags, criterion)
  k <- lo$k
  n <- length(y)
  first <- k + 2L
  if (n - k - 1L < k + 5L) {   # ADF with trend: k + 3 regressors, two residual df
    stop("`lags` = ", k, " is too large for T = ", n,
         " (the regression with constant and trend needs two residual degrees of freedom).",
         call. = FALSE)
  }

  r_c <- adf_regression(y, "c", k, first)
  r_ct <- adf_regression(y, "ct", k, first)
  stats <- c(ADF_c = unname(r_c$tau), ADF_ct = unname(r_ct$tau),
             GLS_c = dfgls_statistic(y, "c", k), GLS_ct = dfgls_statistic(y, "ct", k))
  cv5 <- union_cv5(n, k)
  z <- stats / abs(cv5)
  U <- min(z)
  driver <- names(z)[which.min(z)]

  cv <- rs_critical_values(union_rs, n, k)
  pv <- rs_p_value(union_rs, U, n, k)
  warn <- c(lo$notes, rs_design_warnings(union_rs, n, k))

  other_cv <- lapply(union_components, function(nm) {
    rs <- union_component_rs(nm)
    v <- rs_critical_values(rs, n, k)
    attr(v, "tail") <- "lower"
    attr(v, "p_value") <- rs_p_value(rs, stats[[nm]], n, k)
    attr(v, "source") <- union_component_labels[[nm]]
    v
  })
  names(other_cv) <- union_components

  fit <- r_ct$fit
  coef_table <- cbind(Estimate = fit$coefficients, `Std. Error` = fit$se,
                      `t value` = fit$t_values)
  res <- new_ur_test(
    method = "Union-of-rejections unit-root test (ADF and DF-GLS, constant and trend)",
    test = "Union",
    statistic_name = "U",
    statistic = U,
    null = "unit root",
    alternative = "stationary around a constant or a linear trend",
    trend = "ct",
    lags = k,
    lag_selection = lo$selection,
    nobs = n - k - 1L,
    sample = c(first, n),
    critical_values = cv,
    tail = "lower",
    p_value = pv,
    calibration = list(
      method = "response_surface",
      source = union_source,
      note = paste(c(
        sprintf("finite-sample quantiles of U for T = %d and k = %d lagged differences", n, k),
        sprintf("U = min of t_i / |5%% critical value|, with 5%% critical values %s",
                paste(sprintf("%s %s", union_components, fmt_num(cv5)), collapse = ", ")),
        sprintf("smallest: %s (%s)", driver, union_component_labels[[driver]]),
        p_value_note(pv)), collapse = "; ")),
    level = level,
    other_statistics = stats,
    other_critical_values = other_cv,
    engine = engine_internal(),
    warnings = warn,
    regression = list(coefficients = coef_table, rss = fit$rss,
                      df = fit$df, sigma2 = fit$sigma2,
                      residuals = fit$residuals, X = r_ct$design$X),
    series = series,
    call = call)
  with_residual_warning(res)
}
