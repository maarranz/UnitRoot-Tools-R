# Critical values (spec: "Critical values and calibration")

# ---- tau: MacKinnon response surfaces ---------------------------------------
#
# cv(N) = b_inf + b1 / N + b2 / N^2 + b3 / N^3, evaluated at the effective
# sample size N, for the 1%, 5% and 10% levels (one series, N = 1 in
# MacKinnon's notation). Coefficients for "c" and "ct" are from MacKinnon
# (2010, Table 2); the no-constant case was not updated in 2010 and comes
# from MacKinnon (1996). They agree with statsmodels' adfvalues module.

mackinnon_tau <- list(
  n  = rbind(`0.01` = c(-2.56574, -2.2358,   -3.627,    0),
             `0.05` = c(-1.94100, -0.2686,   -3.365,   31.223),
             `0.1`  = c(-1.61682,  0.2656,   -2.714,   25.364)),
  c  = rbind(`0.01` = c(-3.43035, -6.5393,  -16.786,  -79.433),
             `0.05` = c(-2.86154, -2.8903,   -4.234,  -40.040),
             `0.1`  = c(-2.56677, -1.5384,   -2.809,    0)),
  ct = rbind(`0.01` = c(-3.95877, -9.0531,  -28.428, -134.155),
             `0.05` = c(-3.41049, -4.3904,   -9.036,  -45.374),
             `0.1`  = c(-3.12705, -2.5856,   -3.925,  -22.380))
)

mackinnon_source <- c(
  n  = "MacKinnon (1996) response surface",
  c  = "MacKinnon (2010) response surface",
  ct = "MacKinnon (2010) response surface"
)

tau_critical_values <- function(trend, nobs) {
  b <- mackinnon_tau[[trend]]
  powers <- c(1, 1 / nobs, 1 / nobs^2, 1 / nobs^3)
  drop(b %*% powers)
}

# ---- phi: Dickey and Fuller (1981) tables -------------------------------------
#
# Tables IV-VI of Dickey and Fuller (1981), as reproduced in Enders' Table B.
# Rows: T = 25, 50, 100, 250, 500, infinity; columns: 1%, 5%, 10%.
# Note: urca 1.3-4 repeats the T = 100 values of phi3 in its T = 250 row at
# the 5% and 10% levels; the published values (8.43, 6.34, 5.39) are used here.

df_phi_sizes <- c(25, 50, 100, 250, 500, Inf)

df_phi <- list(
  phi1 = rbind(c(7.88, 5.18, 4.12), c(7.06, 4.86, 3.94), c(6.70, 4.71, 3.86),
               c(6.52, 4.63, 3.81), c(6.47, 4.61, 3.79), c(6.43, 4.59, 3.78)),
  phi2 = rbind(c(8.21, 5.68, 4.67), c(7.02, 5.13, 4.31), c(6.50, 4.88, 4.16),
               c(6.22, 4.75, 4.07), c(6.15, 4.71, 4.05), c(6.09, 4.68, 4.03)),
  phi3 = rbind(c(10.61, 7.24, 5.91), c(9.31, 6.73, 5.61), c(8.73, 6.49, 5.47),
               c(8.43, 6.34, 5.39), c(8.34, 6.30, 5.36), c(8.27, 6.25, 5.34))
)

# Lookup rule (spec): the largest tabulated sample size not exceeding N,
# capped at 500 (the asymptotic row is never used). Below 25 the T = 25 row is
# used and flagged as an extrapolation.
phi_table_row <- function(nobs) {
  finite <- df_phi_sizes[is.finite(df_phi_sizes)]
  eligible <- finite[finite <= nobs]
  if (length(eligible) == 0L) {
    list(row = 1L, size = 25, extrapolated = TRUE)
  } else {
    size <- max(eligible)
    list(row = match(size, df_phi_sizes), size = size, extrapolated = FALSE)
  }
}

phi_critical_values <- function(name, nobs) {
  r <- phi_table_row(nobs)
  cv <- df_phi[[name]][r$row, ]
  names(cv) <- c("0.01", "0.05", "0.1")
  list(values = cv, size = r$size, extrapolated = r$extrapolated)
}
