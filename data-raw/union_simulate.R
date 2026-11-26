# Finite-sample distribution of the union-of-rejections statistic (Decision
# D22; Harvey, Leybourne and Taylor 2012) under the unit-root null (driftless
# random walk), with union_null_stats() in tau_kernel.cpp. Same (T, k) design
# as tau_simulate.R and phi_simulate.R: T = 20..1000; k = 0..min(12,
# Schwert's kmax), plus k = 13..21 for T >= 150; maximal sample t = k+2..T.
#
# The four statistics (ADF-OLS and DF-GLS, with constant and with trend) are
# normalised by the package's own 5% critical values for (T, k) (tau_rs,
# dfgls_rs, through rs_critical_values(), exactly as union_test() does):
#   U = min_i  stat_i / |cv_i(0.05; T, k)|,
# so each normalised test rejects at 5% below -1, and the union rejects
# below the alpha quantile of U. Stored per (T, k, prob): mean and standard
# error across experiments of the simulated quantile of U, and, as a check,
# the marginal 5% rejection rate of each test.
#
# Run from data-raw/ after make_sysdata.R:  Rscript union_simulate.R
# (about 65 minutes on 2 cores).
Rcpp::sourceCpp("tau_kernel.cpp")
pkgload::load_all("..", quiet = TRUE)
ns <- asNamespace("unitroottools")
RNGkind("L'Ecuyer-CMRG")
Tgrid <- c(20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
E <- 25L; R <- 10000L
cores <- max(1L, parallel::detectCores())
set.seed(20261126)
seeds <- vector("list", length(Tgrid)); s <- .Random.seed
for (i in seq_along(Tgrid)) { seeds[[i]] <- s; s <- parallel::nextRNGStream(s) }

union_cv5 <- function(T, k) {
  abs(c(ns$rs_critical_values(ns$tau_rs$c, T, k, levels = 0.05),
        ns$rs_critical_values(ns$tau_rs$ct, T, k, levels = 0.05),
        ns$rs_critical_values(ns$dfgls_rs$c, T, k, levels = 0.05),
        ns$rs_critical_values(ns$dfgls_rs$ct, T, k, levels = 0.05)))
}

one <- function(i) {
  assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
  T <- Tgrid[i]
  lags <- 0:min(12L, as.integer(floor(12 * (T / 100)^0.25)))
  if (T >= 150) lags <- c(lags, 13:21)
  cv <- lapply(lags, function(k) union_cv5(T, k))
  q <- array(NA_real_, c(E, length(lags), length(probs)))
  rej <- array(NA_real_, c(E, length(lags), 4))
  for (e in seq_len(E)) {
    st <- union_null_stats(as.integer(T), as.integer(lags), R)
    for (j in seq_along(lags)) {
      z <- sweep(st[, j, ], 2, cv[[j]], "/")
      u <- apply(z, 1, min)
      q[e, j, ] <- quantile(u, probs, names = FALSE, type = 8, na.rm = TRUE)
      rej[e, j, ] <- colMeans(z < -1, na.rm = TRUE)
    }
  }
  do.call(rbind, lapply(seq_along(lags), function(j)
    data.frame(T = T, k = lags[j], prob = probs,
               q = colMeans(q[, j, ]), se = apply(q[, j, ], 2, sd) / sqrt(E),
               rej_tau2 = mean(rej[, j, 1]), rej_tau3 = mean(rej[, j, 2]),
               rej_gls_c = mean(rej[, j, 3]), rej_gls_ct = mean(rej[, j, 4]))))
}
t0 <- Sys.time()
res <- parallel::mclapply(order(-Tgrid), one, mc.cores = cores, mc.preschedule = FALSE)
# T = 20 with the largest k leaves no residual degrees of freedom: such a job
# fails and is dropped; the surface is fitted on T >= 25 anyway.
res <- Filter(is.data.frame, res)
out <- do.call(rbind, res)
out <- out[order(out$T, out$k, out$prob), ]
write.csv(out, gzfile("union_quantiles.csv.gz"), row.names = FALSE)
cat("done in", format(Sys.time() - t0), "; rows", nrow(out), "\n")
