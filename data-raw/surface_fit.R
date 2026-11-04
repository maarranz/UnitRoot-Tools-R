# Response-surface fitting shared by the calibration scripts (D11, D15, D16).
#
# `d`: simulated quantiles for one statistic (and trend), columns T, k, prob,
# q (mean over experiments), se (its standard error). For each probability,
# weighted least squares (weights 1/se^2) on T >= Tmin of
#   q(T, k) = theta_inf + a_k / T + b_k / T^2 + c_k / T^3
# (a common asymptotic value, one cubic in 1/T for each lag order k).
# Returns the surface object used by R/calibration.R.
fit_surface <- function(d, Tmin = 25) {
  d <- d[d$T >= Tmin, ]
  probs <- sort(unique(d$prob))
  ks <- sort(unique(d$k))
  intercept <- numeric(length(probs))
  slopes <- array(NA_real_, c(length(ks), length(probs), 3),
                  dimnames = list(k = ks, prob = probs, term = c("1/T", "1/T^2", "1/T^3")))
  diag <- data.frame(prob = probs, ratio = NA_real_, maxabs = NA_real_, rms = NA_real_)
  for (i in seq_along(probs)) {
    dd <- d[abs(d$prob - probs[i]) < 1e-12, ]
    dd$w <- 1 / dd$se^2
    dd$kf <- factor(dd$k, levels = ks)
    if (length(ks) > 1L) {
      m <- lm(q ~ kf:I(1/T) + kf:I(1/T^2) + kf:I(1/T^3), data = dd, weights = w)
      cf <- coef(m)
      for (j in seq_along(ks)) {
        lv <- paste0("kf", ks[j])
        slopes[j, i, ] <- cf[c(paste0(lv, ":I(1/T)"), paste0(lv, ":I(1/T^2)"), paste0(lv, ":I(1/T^3)"))]
      }
    } else {                                  # a single lag order (e.g. k = 0)
      m <- lm(q ~ I(1/T) + I(1/T^2) + I(1/T^3), data = dd, weights = w)
      cf <- coef(m)
      slopes[1, i, ] <- cf[c("I(1/T)", "I(1/T^2)", "I(1/T^3)")]
    }
    intercept[i] <- cf[["(Intercept)"]]
    r <- m$residuals
    diag[i, -1] <- c(sum(dd$w * r^2) / m$df.residual, max(abs(r)), sqrt(mean(r^2)))
  }
  design <- aggregate(T ~ k, data = d, FUN = min)
  names(design)[2] <- "Tmin"
  list(probs = probs, intercept = intercept, slopes = slopes, ks = ks,
       design = design, fit = diag)
}

report_fit <- function(rs, label, probs = c(0.01, 0.05, 0.10, 0.90, 0.95, 0.99)) {
  f <- rs$fit[vapply(rs$fit$prob, function(p) any(abs(p - probs) < 1e-12), logical(1)), ]
  cat(sprintf("%-8s max |residual| at p in {%s}: %.4f (RMS <= %.4f)\n", label,
              paste(probs, collapse = ", "), max(f$maxabs), max(f$rms)))
}

# Tabulated calibration (D15), for statistics whose small-sample quantiles
# are too irregular for a smooth surface (the Ng-Perron and ERS statistics,
# whose tails explode at small T with many lags): the simulated quantiles
# themselves, by (k, T, prob), interpolated in 1/T by R/calibration.R; above
# the largest simulated T, interpolation towards the limit theta_inf(p). The
# limit does not depend on k; it is estimated by fit_surface() on the
# well-behaved k = 0 cells with T >= Tlim.
table_calibration <- function(d, Tlim = 50) {
  probs <- sort(unique(d$prob))
  ks <- sort(unique(d$k))
  Tgrid <- sort(unique(d$T))
  tab <- array(NA_real_, c(length(ks), length(Tgrid), length(probs)),
               dimnames = list(k = ks, T = Tgrid, prob = probs))
  idx <- cbind(match(d$k, ks), match(d$T, Tgrid), match(round(d$prob, 6), round(probs, 6)))
  tab[idx] <- d$q
  lim <- fit_surface(d[d$k == 0, ], Tmin = Tlim)
  design <- aggregate(T ~ k, data = d, FUN = min)
  names(design)[2] <- "Tmin"
  list(type = "table", probs = probs, ks = ks, Tgrid = Tgrid, table = tab,
       limit = lim$intercept, design = design, limit_fit = lim$fit)
}
