# Writes inst/fixtures/tau_response_surfaces.json: the tau surfaces (D21) for
# other implementations (Python), in the format of dfgls_response_surfaces.json,
# with reference values computed by the package. Run from data-raw/ after
# make_sysdata.R:  Rscript tau_export_json.R
pkgload::load_all("..", quiet = TRUE)
ns <- asNamespace("unitroottools")
rs <- ns$tau_rs
trend_json <- function(r) list(
  probs = r$probs, k = r$ks,
  design = list(k = r$design$k, Tmin = r$design$Tmin),
  intercept = r$intercept,
  slopes = lapply(seq_along(r$ks), function(a) lapply(seq_along(r$probs), function(b) unname(r$slopes[a, b, ]))))
rows <- list()
for (tr in c("n", "c", "ct")) for (T in c(25L, 30L, 47L, 100L, 262L, 800L)) for (k in c(0L, 2L, 5L, 12L, 21L)) {
  nobs <- T - k - 1L
  if (nobs < 10L) next
  for (stat in c(-4.2, -3.1, -2.2, -1.1)) {
    tc <- ns$tau_calibration(tr, stat, nobs, k)
    rows[[length(rows) + 1]] <- list(trend = tr, T = T, k = k, stat = stat, p = tc$p_value,
                                     cv01 = unname(tc$critical_values[1]), cv05 = unname(tc$critical_values[2]),
                                     cv10 = unname(tc$critical_values[3]))
  }
}
out <- list(
  id = "tau_response_surfaces",
  description = paste("Finite-sample quantiles of the ADF t statistics tau1 (no deterministic terms, trend 'n'),",
                      "tau2 (constant, 'c') and tau3 (constant and trend, 'ct') under the unit-root null, as response",
                      "surfaces in T (length of the series) and k (lagged differences). Decision D21 of the specification."),
  model = "q_p(T, k) = intercept[p] + slopes[k][p][1] / T + slopes[k][p][2] / T^2 + slopes[k][p][3] / T^3, with T = N + k + 1 (N the effective number of observations)",
  usage = list(
    critical_value = "evaluate the model at p = 0.01, 0.05, 0.10",
    p_value = paste("F(stat) = pnorm of the linear interpolation of qnorm(probs) against the quantiles at the statistic",
                    "(monotone, exact at the tabulated probabilities); p = F for this lower-tail test. Beyond the quantiles",
                    "do not extrapolate: report the bound 0.001 (stat <= q[0.001]) or 0.999 (stat >= q[0.999])"),
    outside_design = "use k = min(k, max(k)); use T = max(T, Tmin[k]) from design; never extrapolate; warn",
    phillips_perron = "Z_tau uses the k = 0 surface at T = N + 1"),
  statistic = paste("t ratio of y_{t-1} in the OLS regression of diff(y)_t on the deterministic terms (none, constant,",
                    "or constant and t), y_{t-1} and k lagged differences, t = k+2..T"),
  source = list(
    code = "data-raw/ (tau_kernel.cpp, tau_simulate.R, new_response_surfaces.R, tau_validate.R)",
    simulation = paste("driftless random walk, N(0,1) innovations; T = 20..1000, k = 0..min(12, Schwert kmax) and 13..21",
                       "for T >= 150; 40 experiments x 10,000 replications per T, all k on the same paths; fit on T >= 25",
                       "by WLS with weights 1/se^2"),
    validation = paste("out-of-sample size at 5%: 4.7-5.4% (data-raw/tau_validate.R); at k = 0 within 0.008 of",
                       "MacKinnon (2010) critical values")),
  trends = lapply(rs, trend_json),
  reference_values = list(
    note = sprintf("R unitroottools %s: critical values at 1/5/10%% and p-value of `stat`; an implementation reading this file should reproduce them to 1e-12",
                   as.character(utils::packageVersion("unitroottools"))),
    rows = rows))
jsonlite::write_json(out, "../inst/fixtures/tau_response_surfaces.json", auto_unbox = TRUE, digits = I(17), pretty = FALSE)
cat(length(rows), "reference rows written\n")
