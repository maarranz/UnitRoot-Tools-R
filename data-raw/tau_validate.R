# Size check of the tau response surfaces (D21) on fresh simulations, with
# seeds independent of tau_simulate.R. Run from data-raw/ at the package root's
# data-raw folder:  Rscript tau_validate.R
# For each case, sample size T (series length) and lag order k, the empirical
# rejection rates at 1/5/10% using (a) the critical values and (b) p <= level.
# Sample sizes include values between the simulated grid points.
Rcpp::sourceCpp("tau_kernel.cpp")
pkgload::load_all("..", quiet = TRUE)
ns <- asNamespace("unitroottools")
cases <- c("n", "c", "ct")

size_table <- function(T, ks, R, seed) {
  set.seed(seed)
  st <- tau_null_stats(as.integer(T), as.integer(ks), R)
  out <- list()
  for (j in seq_along(ks)) for (cs in seq_along(cases)) {
    s <- st[, j, cs]; s <- s[is.finite(s)]
    nobs <- T - ks[j] - 1L
    tc <- ns$tau_calibration(cases[cs], 0, nobs, ks[j])
    cv <- unname(tc$critical_values)
    p <- vapply(s, function(x) ns$tau_calibration(cases[cs], x, nobs, ks[j])$p_value, numeric(1))
    out[[length(out) + 1]] <- data.frame(
      case = cases[cs], T = T, k = ks[j],
      cv01 = mean(s < cv[1]), cv05 = mean(s < cv[2]), cv10 = mean(s < cv[3]),
      p01 = mean(p <= 0.01), p05 = mean(p <= 0.05), p10 = mean(p <= 0.10))
  }
  do.call(rbind, out)
}

R <- 20000
grid <- list(list(T = 33L, ks = c(0L, 2L)), list(T = 60L, ks = c(0L, 3L, 6L)),
             list(T = 137L, ks = c(0L, 4L, 12L)), list(T = 290L, ks = c(0L, 8L, 16L)),
             list(T = 610L, ks = c(0L, 20L)))
res <- do.call(rbind, parallel::mclapply(seq_along(grid), function(i)
  size_table(grid[[i]]$T, grid[[i]]$ks, R, 4242L + i), mc.cores = parallel::detectCores()))
cat(sprintf("Empirical size, R = %d replications per cell (MC se at 5%%: %.4f)\n", R, sqrt(.05 * .95 / R)))
print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("\nLargest |size - nominal| at 5%%: critical values %.4f, p-values %.4f\n",
            max(abs(res$cv05 - 0.05)), max(abs(res$p05 - 0.05))))

# Comparison with MacKinnon (2010) at k = 0
cat("\nOwn surfaces vs MacKinnon (2010), k = 0 (1%, 5%, 10%)\n")
for (cs in cases) for (N in c(24L, 49L, 99L, 249L, 499L)) {
  own <- unname(ns$tau_calibration(cs, 0, N, 0L)$critical_values)
  mk <- unname(ns$mackinnon_critical_values(cs, N))
  cat(sprintf("%-2s N = %3d  own %s   MacKinnon %s   max diff %.4f\n", cs, N,
              paste(sprintf("%7.4f", own), collapse = ""), paste(sprintf("%7.4f", mk), collapse = ""),
              max(abs(own - mk))))
}
