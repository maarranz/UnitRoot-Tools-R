# DF-GLS response surfaces (spec, Decision D11).
#
# Inputs: simulated quantiles of the DF-GLS statistic under the null
# (data-raw/dfgls_quantiles.csv.gz, k = 0..12; dfgls_quantiles_ext.csv.gz,
# k = 13..21), produced by dfgls_simulate.R with the kernel dfgls_kernel.cpp.
# Each row: trend, T (length of the series), k (lagged differences), prob,
# q (mean over experiments of the simulated quantile), se (its standard
# error across experiments).
#
# Model, for each trend and probability, with a common asymptotic value and
# one cubic in 1/T for each lag order k:
#   q(T, k) = theta_inf + a_k / T + b_k / T^2 + c_k / T^3,
# estimated by weighted least squares (weights 1/se^2) on T >= 25. The
# paired pattern in k (orders 2j and 2j+1 have nearly equal quantiles in
# small samples) rules out smooth polynomials in k.
#
# Output: data-raw/dfgls_rs.rds (`dfgls_rs`); R/sysdata.rda is assembled by
# make_sysdata.R.

read_q <- function(f) if (file.exists(f)) read.csv(f) else NULL
d <- rbind(read_q("dfgls_quantiles.csv.gz"), read_q("dfgls_quantiles_ext.csv.gz"))
Tmin <- 25
d <- d[d$T >= Tmin, ]
probs <- sort(unique(d$prob))
ks <- sort(unique(d$k))

fit_trend <- function(tr) {
  intercept <- numeric(length(probs))
  slopes <- array(NA_real_, c(length(ks), length(probs), 3),
                  dimnames = list(k = ks, prob = probs, term = c("1/T", "1/T^2", "1/T^3")))
  diag <- data.frame(prob = probs, ratio = NA_real_, maxabs = NA_real_, rms = NA_real_)
  for (i in seq_along(probs)) {
    dd <- d[d$trend == tr & abs(d$prob - probs[i]) < 1e-12, ]
    dd$w <- 1 / dd$se^2
    dd$kf <- factor(dd$k, levels = ks)
    m <- lm(q ~ kf:I(1/T) + kf:I(1/T^2) + kf:I(1/T^3), data = dd, weights = w)
    cf <- coef(m)
    intercept[i] <- cf[["(Intercept)"]]
    for (j in seq_along(ks)) {
      lv <- paste0("kf", ks[j])
      slopes[j, i, ] <- cf[c(paste0(lv, ":I(1/T)"), paste0(lv, ":I(1/T^2)"), paste0(lv, ":I(1/T^3)"))]
    }
    r <- m$residuals
    diag[i, -1] <- c(sum(dd$w * r^2) / m$df.residual, max(abs(r)), sqrt(mean(r^2)))
  }
  design <- aggregate(T ~ k, data = d[d$trend == tr, ], FUN = min)
  names(design)[2] <- "Tmin"
  list(probs = probs, intercept = intercept, slopes = slopes, ks = ks,
       design = design, fit = diag)
}

dfgls_rs <- list(c = fit_trend("c"), ct = fit_trend("ct"))
for (tr in names(dfgls_rs)) {
  f <- dfgls_rs[[tr]]$fit
  cat(sprintf("%-2s: max |residual| %.4f, max RMS %.4f, WRSS/df median %.2f, max %.2f\n",
              tr, max(f$maxabs), max(f$rms), median(f$ratio), max(f$ratio)))
}
if (!interactive()) {
  out <- commandArgs(trailingOnly = TRUE)
  if (length(out)) saveRDS(dfgls_rs, out[1])   # then run make_sysdata.R
}
