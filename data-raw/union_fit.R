# Response surface of the union-of-rejections statistic U (D22), from
# union_simulate.R, fitted by fit_surface() (surface_fit.R) on T >= 25, as the
# other surfaces. Writes union_rs.rds; then run make_sysdata.R.
# U is normalised by the tau and DF-GLS surfaces: if tau_rs or dfgls_rs is
# ever refitted, union_simulate.R must be run again.
source("surface_fit.R")
u <- read.csv("union_quantiles.csv.gz")
u <- u[!is.na(u$q), ]
union_rs <- fit_surface(u[, c("T", "k", "prob", "q", "se")])
report_fit(union_rs, "union", probs = c(0.01, 0.05, 0.10))
saveRDS(union_rs, "union_rs.rds")
