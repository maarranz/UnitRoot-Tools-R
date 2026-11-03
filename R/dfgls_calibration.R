# DF-GLS critical values and p-values (spec: Decision D11).
#
# Own finite-sample response surfaces, estimated by
# data-raw/dfgls_response_surfaces.R from the simulation in
# data-raw/dfgls_simulate.R and stored in R/sysdata.rda as `dfgls_rs`. For
# each trend and tabulated probability p, the p-quantile of the statistic
# with T observations and k lagged differences is
#   q_p(T, k) = theta_inf(p) + a_k(p) / T + b_k(p) / T^2 + c_k(p) / T^3.

# Fitted quantiles at all tabulated probabilities. Outside the simulation
# design the surface is not extrapolated: k is capped at the largest simulated
# order and T is raised to the smallest sample size simulated for that k
# (dfgls_design_warnings() says so).
dfgls_quantiles <- function(trend, T, k) {
  rs <- dfgls_rs[[trend]]
  kk <- min(k, max(rs$ks))
  T <- max(T, rs$design$Tmin[match(kk, rs$design$k)])
  drop(rs$intercept + rs$slopes[match(kk, rs$ks), , ] %*% c(1 / T, 1 / T^2, 1 / T^3))
}

dfgls_critical_values <- function(trend, T, k) {
  rs <- dfgls_rs[[trend]]
  q <- dfgls_quantiles(trend, T, k)
  cv <- vapply(c(0.01, 0.05, 0.10), function(p) q[[which.min(abs(rs$probs - p))]],
               numeric(1))
  names(cv) <- c("0.01", "0.05", "0.1")
  cv
}

# p-value as in MacKinnon (1996): local cubic regression of qnorm(p) on the
# fitted quantiles nearest to the statistic, evaluated at the statistic.
dfgls_p_value <- function(stat, trend, T, k, npoints = 11L) {
  rs <- dfgls_rs[[trend]]
  q <- dfgls_quantiles(trend, T, k)
  idx <- order(abs(q - stat))[seq_len(npoints)]
  qq <- q[idx]
  fit <- stats::lm.fit(cbind(1, qq, qq^2, qq^3), stats::qnorm(rs$probs[idx]))
  p <- stats::pnorm(sum(fit$coefficients * c(1, stat, stat^2, stat^3)))
  min(max(p, 0), 1)
}

# Warnings when (T, k) lies outside the simulation design.
dfgls_design_warnings <- function(trend, T, k) {
  rs <- dfgls_rs[[trend]]
  w <- character()
  if (k > max(rs$ks)) {
    w <- c(w, sprintf(paste("k = %d exceeds the largest simulated lag order (%d);",
                            "critical values and p-value use k = %d."), k, max(rs$ks), max(rs$ks)))
  }
  tmin <- rs$design$Tmin[match(min(k, max(rs$ks)), rs$design$k)]
  if (T < tmin) {
    w <- c(w, sprintf(paste("T = %d is below the smallest sample size simulated for",
                            "k = %d (T = %d); critical values and p-value are those for T = %d."),
                      T, min(k, max(rs$ks)), tmin, tmin))
  }
  w
}
