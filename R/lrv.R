# Long-run variance (spec: Decision D13)
#
# Kernel estimators of the long-run variance of a (mean-zero) series u,
#   lambda^2 = gamma_0 + 2 sum_{j >= 1} w(j / b) gamma_j,
#   gamma_j = n^{-1} sum_{t = j+1}^{n} u_t u_{t-j},
# with the Bartlett, Parzen and Quadratic Spectral kernels (Andrews, 1991),
# and bandwidths fixed, by Schwert's rules, or chosen automatically by
# Andrews (1991, AR(1) plug-in) or Newey and West (1994, nonparametric).

lrv_kernels <- c("bartlett", "parzen", "qs")

kernel_weight <- function(x, kernel) {
  x <- abs(x)
  switch(kernel,
    bartlett = ifelse(x <= 1, 1 - x, 0),
    parzen = ifelse(x <= 0.5, 1 - 6 * x^2 + 6 * x^3,
                    ifelse(x <= 1, 2 * (1 - x)^3, 0)),
    qs = {
      z <- 6 * pi * x / 5
      w <- 25 / (12 * pi^2 * x^2) * (sin(z) / z - cos(z))
      w[x == 0] <- 1
      w
    })
}

autocovariances <- function(u, maxlag) {
  n <- length(u)
  vapply(0:maxlag, function(j) sum(u[(j + 1):n] * u[1:(n - j)]) / n, numeric(1))
}

# Andrews (1991), univariate AR(1) plug-in: alpha(1) for Bartlett, alpha(2)
# for Parzen and QS (eqs. 6.2-6.4).
andrews_bandwidth <- function(u, kernel) {
  n <- length(u)
  rho <- sum(u[-1] * u[-n]) / sum(u[-n]^2)
  rho <- max(min(rho, 0.97), -0.97)      # keep the plug-in finite
  if (kernel == "bartlett") {
    a1 <- 4 * rho^2 / ((1 - rho)^2 * (1 + rho)^2)
    1.1447 * (a1 * n)^(1 / 3)
  } else {
    a2 <- 4 * rho^2 / (1 - rho)^4
    c(parzen = 2.6614, qs = 1.3221)[[kernel]] * (a2 * n)^(1 / 5)
  }
}

# Newey and West (1994), nonparametric automatic bandwidth.
newey_west_bandwidth <- function(u, kernel) {
  n <- length(u)
  m <- switch(kernel,
              bartlett = floor(4 * (n / 100)^(2 / 9)),
              parzen = floor(4 * (n / 100)^(4 / 25)),
              qs = floor(4 * (n / 100)^(2 / 25)))
  m <- min(m, n - 1L)
  g <- autocovariances(u, m)
  j <- seq_len(m)
  s0 <- g[1] + 2 * sum(g[-1])
  if (kernel == "bartlett") {
    s1 <- 2 * sum(j * g[-1])
    1.1447 * ((s1 / s0)^2)^(1 / 3) * n^(1 / 3)
  } else {
    s2 <- 2 * sum(j^2 * g[-1])
    c(parzen = 2.6614, qs = 1.3221)[[kernel]] * ((s2 / s0)^2)^(1 / 5) * n^(1 / 5)
  }
}

# Resolve a `bandwidth` argument. Schwert's rules and a fixed number give a
# number of autocovariances l, used with the Bartlett weights 1 - j/(l+1)
# (Newey-West's convention, as in Kwiatkowski et al. 1992 and urca), i.e.
# b = l + 1; the automatic rules give b directly.
resolve_bandwidth <- function(bandwidth, u, kernel, n_rule = length(u), n_label = "n") {
  if (is.character(bandwidth) && length(bandwidth) == 1L) {
    if (bandwidth %in% c("short", "long")) {
      m <- if (bandwidth == "short") 4 else 12
      l <- as.integer(trunc(m * (n_rule / 100)^(1 / 4)))
      return(list(b = l + 1, lags = l, automatic = FALSE,
                  rule = sprintf("%s rule, trunc(%d (%s/100)^(1/4))", bandwidth, m, n_label)))
    }
    if (bandwidth %in% c("andrews", "newey-west")) {
      b <- if (bandwidth == "andrews") andrews_bandwidth(u, kernel) else newey_west_bandwidth(u, kernel)
      return(list(b = b, lags = NA_integer_, automatic = TRUE,
                  rule = if (bandwidth == "andrews") "Andrews (1991), AR(1) plug-in"
                         else "Newey-West (1994)"))
    }
  }
  if (is.numeric(bandwidth) && !is.logical(bandwidth)) {
    l <- check_order(bandwidth, "bandwidth")
    return(list(b = l + 1, lags = l, automatic = FALSE, rule = "fixed by the user"))
  }
  stop("`bandwidth` must be \"short\", \"long\", \"andrews\", \"newey-west\" or a ",
       "non-negative integer.", call. = FALSE)
}

# Text for the result: what was used, in the user's terms.
bandwidth_text <- function(kernel, bw) {
  kname <- c(bartlett = "Bartlett", parzen = "Parzen", qs = "Quadratic Spectral")[[kernel]]
  if (bw$automatic) {
    sprintf("%s kernel, bandwidth %.3f (%s)", kname, bw$b, bw$rule)
  } else {
    sprintf("%s kernel, %d autocovariance%s (%s)", kname, bw$lags,
            if (bw$lags == 1L) "" else "s", bw$rule)
  }
}

# Long-run variance of u with kernel and bandwidth b (in the kernel's scale:
# weights w(j / b)). Returns gamma_0, the estimate, b and the number of
# autocovariances with non-zero weight.
long_run_variance <- function(u, kernel = "bartlett", b) {
  n <- length(u)
  maxlag <- if (kernel == "qs") n - 1L else min(n - 1L, max(0L, ceiling(b) - 1L))
  g <- autocovariances(u, maxlag)
  w <- if (maxlag > 0) kernel_weight(seq_len(maxlag) / b, kernel) else numeric()
  lrv <- g[1] + 2 * sum(w * g[-1])
  list(gamma0 = g[1], lrv = lrv, b = b, lags_used = as.integer(sum(w != 0)))
}

check_kernel <- function(kernel) {
  if (!is.character(kernel) || length(kernel) != 1L || !kernel %in% lrv_kernels) {
    stop("`kernel` must be one of \"bartlett\", \"parzen\", \"qs\".", call. = FALSE)
  }
  kernel
}
