# Exact limiting distributions of the KPSS statistics (spec: Decision D14).
#
# eta = sum_k lambda_k Z_k^2, Z_k iid N(0,1), with lambda_k the eigenvalues
# of the covariance kernel of the first-level Brownian bridge (eta_mu,
# K(s,t) = min(s,t) - st) or of the second-level bridge (eta_tau,
# K(s,t) = min(s,t) - st - 3st(1-s)(1-t)). Eigenvalues by the Nystrom method
# (midpoint rule, N = 2000 points); for eta_mu they equal 1/(k pi)^2 to six
# significant digits, and their sums equal the means 1/6 and 1/15.
# Probabilities by Imhof's (1961) inversion formula (kpss_upper_prob() in
# R/kpss.R). Writes data-raw/kpss_limit.rds; R/sysdata.rda is assembled by
# make_sysdata.R.
kern <- list(mu  = function(s, t) pmin(s, t) - s * t,
             tau = function(s, t) pmin(s, t) - s * t - 3 * s * t * (1 - s) * (1 - t))
N <- 2000
x <- (seq_len(N) - 0.5) / N
kpss_limit <- lapply(kern, function(K) {
  lam <- eigen(outer(x, x, K) / N, symmetric = TRUE, only.values = TRUE)$values
  sort(lam[lam > 1e-12], decreasing = TRUE)
})
saveRDS(kpss_limit, "kpss_limit.rds")
cat("means:", sum(kpss_limit$mu), sum(kpss_limit$tau), "\n")
