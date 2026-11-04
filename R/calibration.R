# Critical values and p-values from the package's own response surfaces
# (spec: Decisions D11, D15, D16).
#
# A surface `rs` (estimated in data-raw/ and stored in R/sysdata.rda) gives,
# for each tabulated probability p, the p-quantile of a statistic computed
# from a series of length T with k lagged differences:
#   q_p(T, k) = theta_inf(p) + a_k(p) / T + b_k(p) / T^2 + c_k(p) / T^3.
# Fields: probs, intercept, slopes [k, p, 3], ks, design (k, Tmin).
# Objects: dfgls_rs (D11: DF-GLS t, by trend); phi_rs and rho_rs (D16).
#
# A tabulated calibration (type "table": gls_rs, D15) holds the simulated
# quantiles themselves, table[k, T, p] on the grid Tgrid, and their limit
# theta_inf(p); quantiles at other T are interpolated linearly in 1/T
# between neighbouring grid points (beyond the largest, towards the limit).

# Fitted quantiles at all tabulated probabilities. Outside the simulation
# design the surface is not extrapolated: k is capped at the largest
# simulated order and T is raised to the smallest sample size simulated for
# that k (rs_design_warnings() says so).
rs_quantiles <- function(rs, T, k) {
  kk <- min(k, max(rs$ks))
  T <- max(T, rs$design$Tmin[match(kk, rs$design$k)])
  if (identical(rs$type, "table")) return(table_quantiles(rs, T, kk))
  drop(rs$intercept + rs$slopes[match(kk, rs$ks), , ] %*% c(1 / T, 1 / T^2, 1 / T^3))
}

table_quantiles <- function(rs, T, k) {
  tab <- rs$table[match(k, rs$ks), , ]
  ok <- !is.na(tab[, 1])
  Ts <- rs$Tgrid[ok]
  tab <- tab[ok, , drop = FALSE]
  i <- findInterval(T, Ts)                # Ts[i] <= T < Ts[i + 1]
  if (i >= length(Ts)) {                  # beyond the grid: towards the limit
    lo <- tab[length(Ts), ]; hi <- rs$limit
    x0 <- 1 / Ts[length(Ts)]; x1 <- 0
  } else {
    lo <- tab[i, ]; hi <- tab[i + 1L, ]
    x0 <- 1 / Ts[i]; x1 <- 1 / Ts[i + 1L]
  }
  w <- (1 / T - x0) / (x1 - x0)
  unname((1 - w) * lo + w * hi)
}

# Critical values for a lower-tail test (reject for small values) or an
# upper-tail test (reject for large values), at the levels given.
rs_critical_values <- function(rs, T, k, tail = "lower", levels = c(0.01, 0.05, 0.10)) {
  q <- rs_quantiles(rs, T, k)
  probs <- if (tail == "lower") levels else 1 - levels
  cv <- vapply(probs, function(p) q[[which.min(abs(rs$probs - p))]], numeric(1))
  names(cv) <- format_levels(levels)
  cv
}

# p-value as in MacKinnon (1996): a local cubic regression of qnorm(p) on the
# fitted quantiles nearest to the statistic gives the distribution function
# at the statistic; the p-value is F (lower tail) or 1 - F (upper tail).
rs_p_value <- function(rs, stat, T, k, tail = "lower", npoints = 11L) {
  q <- rs_quantiles(rs, T, k)
  idx <- order(abs(q - stat))[seq_len(npoints)]
  qq <- q[idx]
  fit <- stats::lm.fit(cbind(1, qq, qq^2, qq^3), stats::qnorm(rs$probs[idx]))
  F <- stats::pnorm(sum(fit$coefficients * c(1, stat, stat^2, stat^3)))
  p <- if (tail == "lower") F else 1 - F
  min(max(p, 0), 1)
}

# Warnings when (T, k) lies outside the simulation design.
rs_design_warnings <- function(rs, T, k) {
  w <- character()
  if (k > max(rs$ks)) {
    w <- c(w, sprintf(paste("k = %d exceeds the largest simulated lag order (%d);",
                            "critical values and p-value use k = %d."), k, max(rs$ks), max(rs$ks)))
  }
  kk <- min(k, max(rs$ks))
  tmin <- rs$design$Tmin[match(kk, rs$design$k)]
  if (T < tmin) {
    w <- c(w, sprintf(paste("T = %d is below the smallest sample size simulated for",
                            "k = %d (T = %d); critical values and p-value are those for T = %d."),
                      T, kk, tmin, tmin))
  }
  w
}

format_levels <- function(levels) {
  vapply(levels, function(l) as.character(l), character(1))
}

# DF-GLS (D11)
dfgls_quantiles <- function(trend, T, k) rs_quantiles(dfgls_rs[[trend]], T, k)
dfgls_critical_values <- function(trend, T, k) rs_critical_values(dfgls_rs[[trend]], T, k)
dfgls_p_value <- function(stat, trend, T, k, npoints = 11L) {
  rs_p_value(dfgls_rs[[trend]], stat, T, k, npoints = npoints)
}
dfgls_design_warnings <- function(trend, T, k) rs_design_warnings(dfgls_rs[[trend]], T, k)
