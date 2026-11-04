# Response surfaces for the phi statistics and the Dickey-Fuller normalized
# bias used for Z_alpha (D16). Inputs from phi_simulate.R and rho_simulate.R;
# fit by fit_surface() (surface_fit.R) on T >= 25. Writes phi_rs.rds and
# rho_rs.rds; then run make_sysdata.R.
source("surface_fit.R")

phi <- read.csv("phi_quantiles.csv.gz")
phi_rs <- lapply(c(phi1 = "phi1", phi2 = "phi2", phi3 = "phi3"), function(s) {
  rs <- fit_surface(phi[phi$statistic == s, ])
  report_fit(rs, s, probs = c(0.90, 0.95, 0.99))
  rs
})
saveRDS(phi_rs, "phi_rs.rds")

rho <- read.csv("rho_quantiles.csv.gz")
rho_rs <- lapply(c(c = "c", ct = "ct"), function(tr) {
  rs <- fit_surface(rho[rho$trend == tr, ])
  report_fit(rs, paste0("rho_", tr), probs = c(0.01, 0.05, 0.10))
  rs
})
saveRDS(rho_rs, "rho_rs.rds")

# GLS-based statistics (D15): Ng-Perron MZa, MZt, MSB, MPT and ERS PT, from
# gls_simulate.R. Their tails explode at small T with many lags (the AR
# spectral estimate when the lag coefficients sum to nearly one), which no
# smooth surface in 1/T follows, so the simulated quantiles are tabulated and
# interpolated (table_calibration(), surface_fit.R).
gls <- read.csv("gls_quantiles.csv.gz")
gls_rs <- lapply(c(MZa = "MZa", MZt = "MZt", MSB = "MSB", MPT = "MPT", PT = "PT"), function(s) {
  lapply(c(c = "c", ct = "ct"), function(tr) {
    table_calibration(gls[gls$statistic == s & gls$trend == tr, ])
  })
})
saveRDS(gls_rs, "gls_rs.rds")
