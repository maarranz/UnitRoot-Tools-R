# Checks of the union-of-rejections calibration (D22). Run from data-raw/
# after make_sysdata.R:  Rscript union_validate.R   (about 10 minutes, 2 cores)
#
# 1. Out-of-sample size: fresh seeds, sample sizes between the grid points,
#    20,000 replications per cell; rejection at 1/5/10% by the critical value
#    and by p <= level (they must agree).
# 2. Power at T = 200 (500 replications): the union against each component
#    test and the original strategy, for stationary AR(1) series with and
#    without a trend and with a small or a large initial condition.
Rcpp::sourceCpp("tau_kernel.cpp")
pkgload::load_all("..", quiet = TRUE)
ns <- asNamespace("unitroottools")
cores <- parallel::detectCores()

size_cell <- function(T, k, R, seed) {
  set.seed(seed)
  st <- union_null_stats(as.integer(T), as.integer(k), R)[, 1, ]
  u <- apply(sweep(st, 2, abs(ns$union_cv5(T, k)), "/"), 1, min)
  u <- u[is.finite(u)]
  cv <- unname(ns$rs_critical_values(ns$union_rs, T, k))
  p <- vapply(u, function(s) ns$rs_p_value(ns$union_rs, s, T, k), numeric(1))
  data.frame(T = T, k = k, cv01 = mean(u < cv[1]), cv05 = mean(u < cv[2]), cv10 = mean(u < cv[3]),
             p01 = mean(p <= 0.01), p05 = mean(p <= 0.05), p10 = mean(p <= 0.10))
}
cells <- expand.grid(k = c(0L, 3L, 8L), T = c(33L, 60L, 137L, 290L, 610L))
cells <- cells[!(cells$T == 33L & cells$k == 8L), ]
R <- 20000
size <- do.call(rbind, parallel::mclapply(seq_len(nrow(cells)), function(i)
  size_cell(cells$T[i], cells$k[i], R, 7100L + i), mc.cores = cores))
cat(sprintf("1. Empirical size of the union test, R = %d per cell (MC se at 5%%: %.4f)\n",
            R, sqrt(.05 * .95 / R)))
print(format(size, digits = 3), row.names = FALSE)
cat(sprintf("Largest |size - 0.05|: critical values %.4f, p-values %.4f\n\n",
            max(abs(size$cv05 - 0.05)), max(abs(size$p05 - 0.05))))

# 2. Power, T = 200, k = 1 fixed (to compare the tests, not lag rules)
T <- 200; Rp <- 500
dgps <- list(
  `AR(0.9), level, small initial`  = function() as.numeric(stats::filter(rnorm(T), 0.9, "recursive")),
  `AR(0.9), level, large initial`  = function() { u <- stats::filter(rnorm(T), 0.9, "recursive", init = 0); u <- as.numeric(u) + 0.9^(1:T) * 5 / sqrt(1 - 0.81); u },
  `AR(0.9), trend, small initial`  = function() 0.05 * (1:T) + as.numeric(stats::filter(rnorm(T), 0.9, "recursive")),
  `AR(0.95), level, small initial` = function() as.numeric(stats::filter(rnorm(T), 0.95, "recursive")),
  `random walk (size)`             = function() cumsum(rnorm(T)))
pw <- do.call(rbind, parallel::mclapply(seq_len(Rp * length(dgps)), function(i) {
  d <- names(dgps)[(i - 1) %/% Rp + 1]; set.seed(8000 + i); y <- dgps[[d]]()
  u <- union_test(y, lags = 1)
  s <- suppressWarnings(ur_strategy(y, lags = 1))
  data.frame(dgp = d, union = u$reject,
             ADF_c = u$other_statistics[["ADF_c"]] < u$other_critical_values$ADF_c[["0.05"]],
             ADF_ct = u$other_statistics[["ADF_ct"]] < u$other_critical_values$ADF_ct[["0.05"]],
             GLS_c = u$other_statistics[["GLS_c"]] < u$other_critical_values$GLS_c[["0.05"]],
             GLS_ct = u$other_statistics[["GLS_ct"]] < u$other_critical_values$GLS_ct[["0.05"]],
             strategy = isTRUE(s$unit_root_rejected))
}, mc.cores = cores))
cat(sprintf("2. Rejection frequencies at 5%%, T = %d, k = 1, %d replications\n", T, Rp))
agg <- aggregate(cbind(union, ADF_c, ADF_ct, GLS_c, GLS_ct, strategy) ~ dgp, pw, mean)
print(agg[match(names(dgps), agg$dgp), ], row.names = FALSE, digits = 3)
