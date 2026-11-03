# Limiting distributions of the KPSS statistics (eta_mu, eta_tau): T = 2000,
# iid N(0,1), 40 experiments x 25,000 replications. Quantiles at 117
# probabilities with standard errors across experiments.
Rcpp::sourceCpp("kpss_kernel.cpp")   # run from data-raw/; checks kpss_limit.R
set.seed(20261106)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
E <- 40L; R <- 25000L; T <- 2000L
q <- array(NA_real_, c(E, 2, length(probs)))
for (e in seq_len(E)) { s <- kpss_limit_stats(T, R); for (j in 1:2) q[e, j, ] <- quantile(s[, j], probs, names = FALSE, type = 8) }
out <- rbind(data.frame(statistic = "eta_mu", prob = probs, q = colMeans(q[, 1, ]), se = apply(q[, 1, ], 2, sd) / sqrt(E)),
             data.frame(statistic = "eta_tau", prob = probs, q = colMeans(q[, 2, ]), se = apply(q[, 2, ], 2, sd) / sqrt(E)))
write.csv(out, "kpss_quantiles.csv", row.names = FALSE)
print(out[abs(out$prob - 0.90) < 1e-9 | abs(out$prob - 0.95) < 1e-9 | abs(out$prob - 0.975) < 1e-9 | abs(out$prob - 0.99) < 1e-9, ], digits = 4)
