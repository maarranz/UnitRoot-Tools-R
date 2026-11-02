# ADF computational core (spec: "The ADF regression")
#
# Internal building blocks used by adf(), select_lags(), the specification
# strategies and, later, the bootstrap and Monte Carlo layers. They work on
# plain numeric vectors and lists (no lm(), no S4) so that they are cheap to
# call many times.
#
# Conventions:
#   * observations are indexed t = 1, ..., T in the original sample;
#   * the regression with k lagged differences uses t = first, ..., T, with
#     first >= k + 2 (first = k + 2 is the maximal sample; first = kmax + 2
#     is the common sample used for lag selection);
#   * the trend regressor is t itself, so its coefficient has the same
#     meaning in every regression of a procedure.

# Design matrix of the ADF regression
#   Delta y_t = a + gamma t + b y_{t-1} + sum_j beta_j Delta y_{t-j} + e_t
# Columns: const, trend, y_lag1, dy_lag1, ..., dy_lagk (deterministics only
# when included by `trend`).
adf_design <- function(y, trend, k, first = k + 2L) {
  n <- length(y)
  k <- as.integer(k)
  first <- as.integer(first)
  if (is.na(k) || k < 0L) stop("`k` must be a non-negative integer.", call. = FALSE)
  if (first < k + 2L) {
    stop("The sample must start at t >= k + 2 = ", k + 2L, ".", call. = FALSE)
  }
  if (first > n) stop("The sample is empty.", call. = FALSE)

  t <- first:n
  dy <- c(NA_real_, diff(y))               # dy[t] = y[t] - y[t-1]

  cols <- list()
  if (trend %in% c("c", "ct")) cols$const <- rep(1, length(t))
  if (trend == "ct") cols$trend <- as.numeric(t)
  cols$y_lag1 <- y[t - 1L]
  for (j in seq_len(k)) cols[[paste0("dy_lag", j)]] <- dy[t - j]

  X <- do.call(cbind, cols)
  list(response = dy[t], X = X, t = t, trend = trend, k = k)
}

# Least-squares fit by QR. Returns coefficients, standard errors, t ratios,
# RSS, residuals and degrees of freedom. A zero-column X (all coefficients
# restricted to zero) is allowed: then RSS = sum(response^2).
qr_fit <- function(X, response) {
  n <- length(response)
  p <- NCOL(X)
  if (p == 0L) {
    return(list(coefficients = numeric(), se = numeric(), t_values = numeric(),
                rss = sum(response^2), residuals = response,
                nobs = n, df = n, sigma2 = sum(response^2) / n))
  }
  if (n <= p) {
    stop("Not enough observations: N = ", n, " with ", p, " regressors.",
         call. = FALSE)
  }
  q <- qr(X)
  if (q$rank < p) {
    stop("The ADF regressors are perfectly collinear (rank ", q$rank,
         " < ", p, "). Is the series constant or too short?", call. = FALSE)
  }
  coef <- qr.coef(q, response)
  resid <- qr.resid(q, response)
  rss <- sum(resid^2)
  df <- n - p
  sigma2 <- rss / df
  # (X'X)^{-1} from R, undoing any column pivoting
  piv <- q$pivot
  XtXinv <- chol2inv(qr.R(q))
  XtXinv[piv, piv] <- XtXinv
  se <- sqrt(diag(XtXinv) * sigma2)
  names(coef) <- names(se) <- colnames(X)
  list(coefficients = coef, se = se, t_values = coef / se,
       rss = rss, residuals = resid, nobs = n, df = df, sigma2 = sigma2)
}

# Restricted designs for the phi statistics, on the same sample:
#   phi1 (trend "c"):  a = b = 0          -> keep lagged differences only
#   phi2 (trend "ct"): a = gamma = b = 0  -> keep lagged differences only
#   phi3 (trend "ct"): gamma = b = 0      -> keep const and lagged differences
phi_restrictions <- list(
  c  = list(phi1 = c("const", "y_lag1")),
  ct = list(phi2 = c("const", "trend", "y_lag1"),
            phi3 = c("trend", "y_lag1"))
)

f_statistic <- function(rss_r, rss_u, r, n, p_u) {
  ((rss_r - rss_u) / r) / (rss_u / (n - p_u))
}

# ADF regression with tau and phi statistics.
# Returns the unrestricted fit, the design, tau (named tau1/tau2/tau3) and the
# phi statistics that belong to the deterministic specification.
adf_regression <- function(y, trend, k, first = k + 2L) {
  d <- adf_design(y, trend, k, first)
  fit <- qr_fit(d$X, d$response)

  tau_name <- c(n = "tau1", c = "tau2", ct = "tau3")[[trend]]
  tau <- stats::setNames(fit$t_values[["y_lag1"]], tau_name)

  phi <- numeric()
  for (nm in names(phi_restrictions[[trend]])) {
    drop <- phi_restrictions[[trend]][[nm]]
    keep <- setdiff(colnames(d$X), drop)
    fit_r <- qr_fit(d$X[, keep, drop = FALSE], d$response)
    phi[[nm]] <- f_statistic(fit_r$rss, fit$rss, length(drop),
                             fit$nobs, NCOL(d$X))
  }

  list(fit = fit, design = d, tau = tau, phi = phi,
       trend = trend, k = d$k, nobs = fit$nobs,
       sample = c(first = d$t[1L], last = d$t[length(d$t)]))
}
