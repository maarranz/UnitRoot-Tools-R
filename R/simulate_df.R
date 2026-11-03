# Monte Carlo of the Dickey-Fuller distribution (teaching tool)

#' Simulate the Dickey-Fuller distribution
#'
#' Generates `R` random walks \eqn{y_t = y_{t-1} + \varepsilon_t},
#' \eqn{y_0 = 0}, \eqn{\varepsilon_t \sim N(0, 1)}, of length `T`, runs the
#' Dickey-Fuller regression (no lagged differences) with the deterministic
#' terms given by `trend`, and records the t ratio of \eqn{b}
#' (\eqn{\tau}) and the normalized coefficient \eqn{N \hat b}.
#'
#' The empirical quantiles of \eqn{\tau} are compared with MacKinnon's
#' response-surface critical values (as used by [adf()]) and with the
#' quantiles of the standard normal distribution: under the unit-root null
#' the t ratio does not have a normal (or Student t) distribution.
#'
#' @param T Length of each simulated series.
#' @param R Number of replications.
#' @param trend `"n"`, `"c"` or `"ct"`: deterministic terms in the test
#'   regression (the data are always a driftless random walk).
#' @param seed Optional seed, for reproducible results. The global random
#'   number generator state is restored on exit.
#'
#' @return An object of class `df_simulation`, with the simulated `tau` and
#'   `normalized_b` values and a `quantiles` table; it has `print` and `plot`
#'   methods.
#' @examples
#' s <- simulate_df(T = 100, R = 2000, trend = "c", seed = 1)
#' s
#' plot(s)
#' @export
simulate_df <- function(T = 100, R = 10000, trend = "c", seed = NULL) {
  trend <- match_trend(trend)
  if (length(T) != 1L || is.na(T) || T < 10 || T != round(T)) {
    stop("`T` must be an integer of at least 10.", call. = FALSE)
  }
  if (length(R) != 1L || is.na(R) || R < 100 || R != round(R)) {
    stop("`R` must be an integer of at least 100.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  tau <- numeric(R)
  nb <- numeric(R)
  for (r in seq_len(R)) {
    y <- cumsum(stats::rnorm(T))
    reg <- adf_regression(y, trend, 0L)
    tau[r] <- reg$tau[[1]]
    nb[r] <- reg$nobs * reg$fit$coefficients[["y_lag1"]]
  }

  nobs <- T - 1L
  levels <- c(0.01, 0.05, 0.10)
  q <- data.frame(
    level = levels,
    simulated = unname(stats::quantile(tau, levels, type = 8)),
    mackinnon = unname(tau_critical_values(trend, nobs)),
    normal = stats::qnorm(levels)
  )
  structure(list(tau = tau, normalized_b = nb, T = T, R = R, trend = trend,
                 nobs = nobs, quantiles = q, seed = seed,
                 statistic_name = c(n = "tau1", c = "tau2", ct = "tau3")[[trend]]),
            class = "df_simulation")
}

#' @export
print.df_simulation <- function(x, ...) {
  cat(sprintf("Simulated Dickey-Fuller distribution of %s\n\n", x$statistic_name))
  cat(sprintf("  Data: driftless random walk, T = %d (N = %d); R = %d replications\n",
              x$T, x$nobs, x$R))
  cat(sprintf("  Test regression: %s\n\n", trend_labels[[x$trend]]))
  q <- x$quantiles
  cat(sprintf("  %-6s %10s %18s %10s\n", "Level", "Simulated", "MacKinnon (resp.)", "N(0,1)"))
  for (i in seq_len(nrow(q))) {
    cat(sprintf("  %-6s %10s %18s %10s\n", fmt_level(q$level[i]),
                fmt_num(q$simulated[i]), fmt_num(q$mackinnon[i]), fmt_num(q$normal[i])))
  }
  se <- sqrt(0.05 * 0.95 / x$R)
  cat(sprintf("\n  Share of simulated tau below the normal 5%% quantile: %.3f",
              mean(x$tau < stats::qnorm(0.05))))
  cat(sprintf("\n  (a correctly sized test would reject 0.050 +/- %.3f)\n", 1.96 * se))
  invisible(x)
}

#' @export
plot.df_simulation <- function(x, ...) {
  rng <- range(c(x$tau, -4.5, 3))
  dens <- stats::density(x$tau)
  h <- graphics::hist(x$tau, breaks = 60, plot = FALSE)
  peak <- max(dens$y, h$density, stats::dnorm(0))
  top <- 1.35 * peak                                     # room for the legend
  graphics::hist(x$tau, breaks = 60, freq = FALSE, xlim = rng, ylim = c(0, top),
                 col = "grey85", border = "white",
                 main = sprintf("%s under the unit-root null", x$statistic_name),
                 sub = sprintf("T = %d, R = %d replications", x$T, x$R),
                 xlab = x$statistic_name, ...)
  graphics::lines(dens, lwd = 2)
  graphics::curve(stats::dnorm(x), add = TRUE, lty = 2, lwd = 2)
  q5 <- x$quantiles$simulated[x$quantiles$level == 0.05]
  graphics::segments(q5, 0, q5, 1.05 * peak, lwd = 2)
  graphics::segments(stats::qnorm(0.05), 0, stats::qnorm(0.05), 1.05 * peak,
                     lty = 2, lwd = 2)
  graphics::legend("topright", bty = "n", lty = c(1, 2), lwd = 2, cex = 0.85,
                   seg.len = 2.5,
                   legend = c(sprintf("simulated, 5%%: %.2f", q5),
                              sprintf("N(0,1), 5%%: %.2f", stats::qnorm(0.05))))
  invisible(x)
}
