# MAIC size/power trade-off quoted in ?select_lags (T = 200, 500 replications,
# 2 cores, ~5 minutes). Run from the package root: Rscript data-raw/maic_tradeoff.R
# Part 1: power (stationary series) and size (unit root) of adf(trend = "c")
# and ur_strategy() with MAIC, AIC, BIC.
pkgload::load_all(quiet=TRUE)
R <- 500; T <- 200
dgps <- list(
  `white noise (diff of RW)` = function() rnorm(T),
  `AR(1) 0.5`  = function() as.numeric(arima.sim(list(ar = 0.5), T)),
  `AR(1) 0.9`  = function() as.numeric(arima.sim(list(ar = 0.9), T)),
  `MA(1) -0.5` = function() as.numeric(arima.sim(list(ma = -0.5), T)),
  `random walk (size)` = function() cumsum(rnorm(T)),
  `ARIMA(1,1,0) 0.5 (size)` = function() cumsum(as.numeric(arima.sim(list(ar = 0.5), T))))
one <- function(d, seed) {
  set.seed(seed); y <- dgps[[d]]()
  out <- NULL
  for (cr in c("MAIC", "AIC", "BIC")) {
    a <- suppressWarnings(adf(y, trend = "c", criterion = cr))
    s <- suppressWarnings(ur_strategy(y, criterion = cr))
    out <- rbind(out, data.frame(dgp = d, crit = cr, k = a$lags, adf_rej = a$reject,
                                 strat_rej = isTRUE(s$unit_root_rejected)))
  }
  out
}
res <- do.call(rbind, parallel::mclapply(seq_len(R * length(dgps)), function(i) {
  d <- names(dgps)[(i - 1) %/% R + 1]; one(d, 1000 + i) }, mc.cores = 2))
agg <- aggregate(cbind(k, adf_rej, strat_rej) ~ dgp + crit, res, mean)
agg$k_median <- aggregate(k ~ dgp + crit, res, median)$k
agg$k_ge12 <- aggregate(k ~ dgp + crit, res, function(v) mean(v >= 12))$k
agg <- agg[order(match(agg$dgp, names(dgps)), agg$crit), ]
print(format(agg, digits = 3), row.names = FALSE)

# Part 2: size under a unit root with negative MA(1) errors.
R <- 500; T <- 200
res <- do.call(rbind, parallel::mclapply(1:(2*R), function(i) {
  th <- if (i <= R) -0.5 else -0.8; set.seed(5000 + i)
  y <- cumsum(as.numeric(arima.sim(list(ma = th), T)))
  do.call(rbind, lapply(c("MAIC","AIC","BIC"), function(cr) {
    a <- suppressWarnings(adf(y, trend = "c", criterion = cr))
    data.frame(theta = th, crit = cr, k = a$lags, rej = a$reject) }))
}, mc.cores = 2))
print(aggregate(cbind(k, rej) ~ theta + crit, res, mean), digits = 3)
