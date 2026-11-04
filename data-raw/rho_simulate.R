# Finite-sample distribution of the Dickey-Fuller normalized-bias statistic
# n * b_hat under the unit-root null (no lagged differences), constant and
# constant-and-trend regressions; used for Phillips-Perron's Z_alpha, which
# has the same limit (Decision D16). T = 20..1000; 40 experiments x 25,000
# replications per T; quantiles at 117 probabilities.
Rcpp::sourceCpp("rho_kernel.cpp")   # run from data-raw/
set.seed(20261108)
Tgrid <- c(20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
E <- 40L; R <- 25000L
out <- NULL
for (T in Tgrid) {
  q <- array(NA_real_, c(E, 2, length(probs)))
  for (e in seq_len(E)) { s <- rho_null_stats(as.integer(T), R); for (j in 1:2) q[e, j, ] <- quantile(s[, j], probs, names = FALSE, type = 8) }
  for (j in 1:2) out <- rbind(out, data.frame(trend = c("c", "ct")[j], T = T, k = 0L, prob = probs,
                                               q = colMeans(q[, j, ]), se = apply(q[, j, ], 2, sd) / sqrt(E)))
}
write.csv(out, gzfile("rho_quantiles.csv.gz"), row.names = FALSE)
cat("done; rows", nrow(out), "\n")
