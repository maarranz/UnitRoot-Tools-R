# Finite-sample distribution of the phi statistics (and tau2, tau3) of the
# ADF regression under the unit-root null (driftless random walk; the phi3
# statistic is invariant to a drift). Design as for DF-GLS (D11): T = 20..1000;
# k = 0..min(12, Schwert's kmax), plus k = 13..21 for T >= 150; maximal sample
# t = k+2..T; E experiments of R replications per T, all k on the same paths.
# Stored per (statistic, T, k, prob): mean and standard error across
# experiments of the simulated quantile.
Rcpp::sourceCpp("phi_kernel.cpp")   # run from data-raw/
RNGkind("L'Ecuyer-CMRG")
Tgrid <- c(20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
stats <- c("phi1", "phi2", "phi3", "tau2", "tau3")
E <- 40L; R <- 10000L
set.seed(20261105)
seeds <- vector("list", length(Tgrid)); s <- .Random.seed
for (i in seq_along(Tgrid)) { seeds[[i]] <- s; s <- parallel::nextRNGStream(s) }
one <- function(i) {
  assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
  T <- Tgrid[i]
  lags <- 0:min(12L, as.integer(floor(12 * (T / 100)^0.25)))
  if (T >= 150) lags <- c(lags, 13:21)
  q <- array(NA_real_, c(E, length(lags), length(stats), length(probs)))
  for (e in seq_len(E)) {
    st <- phi_null_stats(as.integer(T), lags, R)
    for (j in seq_along(lags)) for (s in seq_along(stats))
      q[e, j, s, ] <- quantile(st[, j, s], probs, names = FALSE, type = 8)
  }
  do.call(rbind, lapply(seq_along(stats), function(s) do.call(rbind, lapply(seq_along(lags), function(j)
    data.frame(statistic = stats[s], T = T, k = lags[j], prob = probs,
               q = colMeans(q[, j, s, ]), se = apply(q[, j, s, ], 2, sd) / sqrt(E))))))
}
t0 <- Sys.time()
res <- parallel::mclapply(order(-Tgrid), one, mc.cores = 2, mc.preschedule = FALSE)
# T = 20 with the largest k leaves no residual degrees of freedom (NaN): that
# job fails and is dropped; the surfaces are fitted on T >= 25 anyway.
res <- Filter(is.data.frame, res)
out <- do.call(rbind, res)
out <- out[order(out$statistic, out$T, out$k, out$prob), ]
write.csv(out, gzfile("phi_quantiles.csv.gz"), row.names = FALSE)
cat("done in", format(Sys.time() - t0), "; rows", nrow(out), "\n")
