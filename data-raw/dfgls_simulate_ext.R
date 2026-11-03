# Finite-sample distribution of the DF-GLS statistic under the null.
# Extension: k = 13..21 at every T in Tgrid (main run: dfgls_simulate.R, k <= 12).
# per (trend, T): E experiments of R replications; all lag orders computed on
# the same paths. Stored per cell: mean and standard error (across
# experiments) of the quantiles at probs.
Rcpp::sourceCpp("dfgls_kernel.cpp")   # run from data-raw/
RNGkind("L'Ecuyer-CMRG")
Tgrid <- c(150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
E <- 40L; R <- 10000L
jobs <- expand.grid(T = Tgrid, trend = c("c", "ct"), stringsAsFactors = FALSE)
set.seed(20261104)
seeds <- vector("list", nrow(jobs)); s <- .Random.seed
for (i in seq_len(nrow(jobs))) { seeds[[i]] <- s; s <- parallel::nextRNGStream(s) }
one <- function(i) {
  assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
  T <- jobs$T[i]; tr <- jobs$trend[i]
  lags <- 13:21
  q <- array(NA_real_, c(E, length(lags), length(probs)))
  for (e in seq_len(E)) {
    st <- dfgls_null_stats(as.integer(T), lags, tr == "ct", R)
    for (j in seq_along(lags)) q[e, j, ] <- quantile(st[, j], probs, names = FALSE, type = 8)
  }
  do.call(rbind, lapply(seq_along(lags), function(j) data.frame(
    trend = tr, T = T, k = lags[j], prob = probs,
    q = colMeans(q[, j, ]), se = apply(q[, j, ], 2, sd) / sqrt(E))))
}
t0 <- Sys.time()
res <- parallel::mclapply(order(-jobs$T), one, mc.cores = 2, mc.preschedule = FALSE)
out <- do.call(rbind, res)
out <- out[order(out$trend, out$T, out$k, out$prob), ]
write.csv(out, gzfile("dfgls_quantiles_ext.csv.gz"), row.names = FALSE)
cat("done in", format(Sys.time() - t0), "; rows", nrow(out), "\n")
