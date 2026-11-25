# Finite-sample distribution of the ADF t statistics tau1 (no constant), tau2
# (constant) and tau3 (constant and trend) under the unit-root null (driftless
# random walk), with the kernel in tau_kernel.cpp (nvars = 1). Same design as
# phi_simulate.R (D11, D16): T = 20..1000; k = 0..min(12, Schwert's kmax),
# plus k = 13..21 for T >= 150; maximal sample t = k+2..T; E experiments of R
# replications per T, all k on the same paths. Stored per (statistic, T, k,
# prob): mean and standard error across experiments of the simulated quantile.
# The tau2 and tau3 quantiles of phi_quantiles.csv.gz (independent seeds) are
# a check of this run.
#
# The kernel takes `nvars`: the same script with nvars > 1 gives Engle-Granger
# residual-based statistics (cointegrationtools), cases n, c, ct referring to
# the cointegrating regression.
Rcpp::sourceCpp("tau_kernel.cpp")   # run from data-raw/
RNGkind("L'Ecuyer-CMRG")
nvars <- 1L
Tgrid <- c(20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
stats <- c("tau1", "tau2", "tau3")
E <- 40L; R <- 10000L
cores <- max(1L, parallel::detectCores() - 0L)
set.seed(20261125)
seeds <- vector("list", length(Tgrid)); s <- .Random.seed
for (i in seq_along(Tgrid)) { seeds[[i]] <- s; s <- parallel::nextRNGStream(s) }
one <- function(i) {
  assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
  T <- Tgrid[i]
  lags <- 0:min(12L, as.integer(floor(12 * (T / 100)^0.25)))
  if (T >= 150) lags <- c(lags, 13:21)
  q <- array(NA_real_, c(E, length(lags), length(stats), length(probs)))
  for (e in seq_len(E)) {
    st <- tau_null_stats(as.integer(T), as.integer(lags), R, nvars)
    for (j in seq_along(lags)) for (s in seq_along(stats))
      q[e, j, s, ] <- quantile(st[, j, s], probs, names = FALSE, type = 8, na.rm = TRUE)
  }
  do.call(rbind, lapply(seq_along(stats), function(s) do.call(rbind, lapply(seq_along(lags), function(j)
    data.frame(statistic = stats[s], T = T, k = lags[j], prob = probs,
               q = colMeans(q[, j, s, ]), se = apply(q[, j, s, ], 2, sd) / sqrt(E))))))
}
t0 <- Sys.time()
res <- parallel::mclapply(order(-Tgrid), one, mc.cores = cores, mc.preschedule = FALSE)
# T = 20 with the largest k leaves no residual degrees of freedom: such a job
# fails and is dropped; the surfaces are fitted on T >= 25 anyway.
res <- Filter(is.data.frame, res)
out <- do.call(rbind, res)
out <- out[order(out$statistic, out$T, out$k, out$prob), ]
write.csv(out, gzfile("tau_quantiles.csv.gz"), row.names = FALSE)
cat("done in", format(Sys.time() - t0), "; rows", nrow(out), "\n")
