# Adapter to bootUR (spec: Decision D17): bootstrap unit-root tests by
# Smeekes and Wilms (2023), returned as `ur_test` results.

engine_bootur <- function() {
  paste("bootUR", as.character(utils::packageVersion("bootUR")))
}

bootur_methods <- c(AWB = "autoregressive wild bootstrap",
                    BWB = "block wild bootstrap",
                    DWB = "dependent wild bootstrap",
                    MBB = "moving block bootstrap",
                    SB = "sieve bootstrap",
                    SWB = "sieve wild bootstrap")

#' Bootstrap unit-root tests (bootUR)
#'
#' Bootstrap ADF and union tests computed by the package bootUR (Smeekes and
#' Wilms, 2023), returned in the common result format. bootUR must be
#' installed (`install.packages("bootUR")`); it is not needed for anything
#' else in this package.
#'
#' * `union = FALSE`: bootstrap ADF test with the deterministic terms of
#'   `trend`, on the series detrended by OLS (`detrend = "OLS"`) or by GLS
#'   (`"GLS"`, bootUR's "QD").
#' * `union = TRUE`: the union-of-rejections test of Harvey, Leybourne and
#'   Taylor (2012) as bootstrapped by Smeekes and Taylor (2012): the four
#'   tests with constant or constant and trend, OLS or GLS detrending; `trend`
#'   and `detrend` are not used.
#'
#' The bootstrap gives a p-value, not critical values; the decision is taken
#' from the p-value. bootUR computes its own statistic (detrending first,
#' then a regression without deterministic terms), so it can differ
#' slightly from [adf()] on the same data; the result reports bootUR's.
#'
#' Lags are selected by `criterion` over \eqn{0, \dots, k_{\max}}; by default
#' bootUR rescales the criteria as Cavaliere et al. (2015) propose
#' (`criterion_scale = TRUE`), which is robust to heteroskedasticity. bootUR
#' 1.0.5 cannot fix the number of lags above 0 (it fails when its minimum lag
#' is 1 or more), so `lags` may only be 0 here.
#'
#' Results depend on the random-number stream: give `seed` (or call
#' [set.seed()]) for reproducible p-values.
#'
#' @inheritParams adf
#' @param bootstrap Bootstrap method: `"AWB"` (default), `"BWB"`, `"DWB"`,
#'   `"MBB"`, `"SB"` or `"SWB"` (see bootUR's `boot_adf()`).
#' @param B Number of bootstrap replications (default 1999).
#' @param union Union-of-rejections test instead of a single ADF test.
#' @param detrend `"OLS"` (default) or `"GLS"`.
#' @param criterion_scale Rescaled information criteria (bootUR's default).
#' @param lags `NULL` (select) or 0.
#' @param seed Optional seed for the bootstrap.
#' @param level Significance level of the reported decision.
#' @return An object of class `ur_test`.
#' @references
#' Cavaliere, G., Phillips, P. C. B., Smeekes, S. and Taylor, A. M. R. (2015).
#' Lag length selection for unit root tests in the presence of nonstationary
#' volatility. *Econometric Reviews*, 34, 512-536.
#'
#' Harvey, D. I., Leybourne, S. J. and Taylor, A. M. R. (2012). Testing for
#' unit roots in the presence of uncertainty over both the trend and initial
#' condition. *Journal of Econometrics*, 169, 188-195.
#'
#' Smeekes, S. and Taylor, A. M. R. (2012). Bootstrap union tests for unit
#' roots in the presence of nonstationary volatility. *Econometric Theory*,
#' 28, 422-456.
#'
#' Smeekes, S. and Wilms, I. (2023). bootUR: An R package for bootstrap unit
#' root tests. *Journal of Statistical Software*, 106(12), 1-39.
#' @examples
#' if (requireNamespace("bootUR", quietly = TRUE)) {
#'   y <- cumsum(rnorm(120))
#'   boot_ur(y, trend = "c", B = 199, seed = 1)
#' }
#' @export
boot_ur <- function(y, trend = "c", bootstrap = c("AWB", "BWB", "DWB", "MBB", "SB", "SWB"),
                    B = 1999, union = FALSE, detrend = c("OLS", "GLS"),
                    lags = NULL, max_lags = NULL,
                    criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                    criterion_scale = TRUE, seed = NULL, level = 0.05, var = NULL) {
  call <- match.call()
  if (!requireNamespace("bootUR", quietly = TRUE)) {
    stop("boot_ur() needs the package bootUR: install.packages(\"bootUR\").",
         call. = FALSE)
  }
  bootstrap <- match.arg(bootstrap)
  detrend <- match.arg(detrend)
  criterion <- match.arg(criterion)
  check_level(level)
  check_flag(union, "union")
  check_flag(criterion_scale, "criterion_scale")
  B <- check_order(B, "B")
  if (B < 99L) stop("`B` must be at least 99.", call. = FALSE)
  y <- check_series(y, var = var)
  trend <- match_trend(trend, allowed = c("c", "ct"))
  n <- length(y)
  min_lag <- 0L
  max_lag <- if (is.null(max_lags)) default_kmax(n) else check_order(max_lags, "max_lags")
  if (!is.null(lags)) {
    if (check_order(lags, "lags") != 0L) {
      stop("bootUR 1.0.5 cannot fix the number of lags above 0; leave `lags = NULL` ",
           "to select them (use `max_lags` to bound the search).", call. = FALSE)
    }
    max_lag <- 0L
  }
  if (!is.null(seed)) set.seed(check_order(seed, "seed"))
  args <- list(data = y, data_name = "y", bootstrap = bootstrap, B = B,
               min_lag = min_lag, max_lag = max_lag, criterion = criterion,
               criterion_scale = criterion_scale, show_progress = FALSE,
               do_parallel = FALSE)
  if (union) {
    res <- do.call(bootUR::boot_union, args)
  } else {
    args$deterministics <- c(c = "intercept", ct = "trend")[[trend]]
    args$detrend <- c(OLS = "OLS", GLS = "QD")[[detrend]]
    res <- do.call(bootUR::boot_adf, args)
  }
  sel <- res$details$`selected lags`
  k <- if (union) NA_integer_ else as.integer(sel[1])
  stat <- unname(res$statistic)
  p <- unname(res$p.value)
  others <- if (union) {
    s <- res$details$`individual statistics`
    lab <- colnames(s)
    nm <- paste0("t_", ifelse(grepl("trend", lab), "ct", "c"), "_",
                 ifelse(grepl("^QD", lab), "GLS", "OLS"))
    stats::setNames(as.numeric(s[1, ]), nm)
  } else numeric()
  new_ur_test(
    method = if (union) sprintf("Bootstrap union test (%s)", bootUR_name(bootstrap))
             else sprintf("Bootstrap ADF test (%s, %s detrending)", bootUR_name(bootstrap), detrend),
    test = if (union) "Union" else "Bootstrap ADF",
    statistic_name = if (union) "union" else "t_b",
    statistic = stat, null = "unit root",
    alternative = if (union) "stationary" else adf_alternatives[[trend]],
    trend = if (union) "ct" else trend,
    lags = if (union) NA_integer_ else k,
    lag_selection = if (is.null(lags) && !union) {
      list(method = paste0(criterion, if (criterion_scale) " (rescaled)" else ""),
           kmax = max_lag, detrending = detrend)
    } else NULL,
    nobs = if (union) n else n - k - 1L,
    sample = if (union) c(1L, n) else c(k + 2L, n),
    critical_values = numeric(), tail = "lower", p_value = p,
    calibration = list(method = "bootstrap",
                       source = sprintf("bootUR %s, %s, B = %d", as.character(utils::packageVersion("bootUR")),
                                        bootstrap, B),
                       note = if (is.null(seed)) "p-value depends on the random-number stream (no seed given)"
                              else sprintf("seed = %d", seed)),
    level = level, other_statistics = others, engine = engine_bootur(),
    warnings = if (union) sprintf(paste("union test: trend and detrend are not used; lags",
                                        "selected by %s (kmax = %d) for each test: %s"),
                                  criterion, max_lag, paste(as.integer(sel[1, ]), collapse = ", "))
               else character(),
    call = call)
}

bootUR_name <- function(b) bootur_methods[[b]]
