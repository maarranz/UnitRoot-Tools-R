# Writes inst/fixtures/union_response_surface.json: the surface of the union
# statistic U (D22) for other implementations, in the format of
# dfgls_response_surfaces.json (a single "trend", "union"), with reference
# values computed by the package. U is normalised by the tau surfaces
# (tau_response_surfaces.json, trends c and ct) and the DF-GLS surfaces
# (dfgls_response_surfaces.json, trends c and ct). Run from data-raw/ after
# make_sysdata.R:  Rscript union_export_json.R
pkgload::load_all("..", quiet = TRUE)
ns <- asNamespace("unitroottools")
r <- ns$union_rs
surface <- list(
  probs = r$probs, k = r$ks,
  design = list(k = r$design$k, Tmin = r$design$Tmin),
  intercept = r$intercept,
  slopes = lapply(seq_along(r$ks), function(a) lapply(seq_along(r$probs), function(b) unname(r$slopes[a, b, ]))))
rows <- list()
for (T in c(25L, 30L, 47L, 100L, 262L, 800L)) for (k in c(0L, 2L, 5L, 12L, 21L)) {
  if (T - k - 1L < k + 5L) next
  cv5 <- ns$union_cv5(T, k)
  for (U in c(-2.2, -1.4, -1.1, -0.6)) {
    cv <- unname(ns$rs_critical_values(r, T, k))
    rows[[length(rows) + 1]] <- list(trend = "union", T = T, k = k, stat = U,
                                     p = ns$rs_p_value(r, U, T, k),
                                     cv01 = cv[1], cv05 = cv[2], cv10 = cv[3],
                                     cv5_ADF_c = unname(cv5[["ADF_c"]]), cv5_ADF_ct = unname(cv5[["ADF_ct"]]),
                                     cv5_GLS_c = unname(cv5[["GLS_c"]]), cv5_GLS_ct = unname(cv5[["GLS_ct"]]))
  }
}
out <- list(
  id = "union_response_surface",
  description = paste("Finite-sample quantiles of the union-of-rejections statistic U (Harvey, Leybourne and Taylor",
                      "2012) under the unit-root null, as a response surface in T (length of the series) and k",
                      "(lagged differences, common to the four tests). Decision D22 of the specification."),
  model = "q_p(T, k) = intercept[p] + slopes[k][p][1] / T + slopes[k][p][2] / T^2 + slopes[k][p][3] / T^3",
  statistic = paste("U = min over i of stat_i / |cv_i(0.05; T, k)|, i = ADF with constant (tau, trend c),",
                    "ADF with constant and trend (tau, ct), DF-GLS with constant (c), DF-GLS with constant and trend (ct);",
                    "all with k lagged differences on the maximal sample t = k+2..T. cv_i are the 5% critical values of",
                    "tau_response_surfaces.json (trends c, ct) and dfgls_response_surfaces.json (trends c, ct) at the",
                    "same T and k; the reference rows give them (cv5_*)."),
  usage = list(
    critical_value = "evaluate the model at p = 0.01, 0.05, 0.10; reject the unit root when U is below it",
    p_value = paste("lower tail, as for tau: linear interpolation of qnorm(probs) against the quantiles; exactly the",
                    "tabulated probability at a quantile; bounds 0.001/0.999 beyond them"),
    outside_design = "use k = min(k, max(k)); use T = max(T, Tmin[k]) from design; never extrapolate; warn",
    lags = "one k for the four tests, selected (by default MAIC) on the ADF regression with constant and trend"),
  source = list(
    code = "data-raw/ (tau_kernel.cpp union_null_stats, union_simulate.R, union_fit.R, union_validate.R)",
    simulation = "driftless random walk, N(0,1); same (T, k) design as tau; 25 x 10,000 replications per T"),
  trends = list(union = surface),
  reference_values = list(
    note = sprintf("R unitroottools %s: critical values at 1/5/10%% and p-value of `stat`; to 1e-12",
                   as.character(utils::packageVersion("unitroottools"))),
    rows = rows))
jsonlite::write_json(out, "../inst/fixtures/union_response_surface.json", auto_unbox = TRUE,
                     digits = I(17), pretty = FALSE)
cat(length(rows), "reference rows written\n")
