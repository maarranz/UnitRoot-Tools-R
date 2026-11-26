# p-values equal the level exactly at the tabulated critical values (ported
# from cointegrationtools 0.1.1, audit M03). Before, pnorm(qnorm(p)) rounding
# (and 1 - 0.90 in the upper tail, and uniroot() for KPSS) could give p < level
# at the level's own critical value, contradicting the decision.

test_that("response surfaces: p is exactly the level at the critical value", {
  sets <- list(
    list(tau_rs$n, "lower"), list(tau_rs$c, "lower"), list(tau_rs$ct, "lower"),
    list(dfgls_rs$c, "lower"), list(dfgls_rs$ct, "lower"),
    list(phi_rs$phi1, "upper"), list(phi_rs$phi2, "upper"), list(phi_rs$phi3, "upper"))
  for (s in sets) for (T in c(25L, 77L, 400L)) for (k in c(0L, 3L, 6L)) {
    for (a in c(0.01, 0.05, 0.10)) {
      cv <- unname(rs_critical_values(s[[1]], T, k, tail = s[[2]], levels = a))
      expect_identical(rs_p_value(s[[1]], cv, T, k, tail = s[[2]]), a)
    }
  }
  # still continuous: just inside the critical value, p is just below the level
  cv <- unname(rs_critical_values(tau_rs$c, 100, 2, levels = 0.05))
  p <- rs_p_value(tau_rs$c, cv - 1e-6, 100, 2)
  expect_lt(p, 0.05)
  expect_gt(p, 0.0499)
})

test_that("KPSS: p is exactly the level at the critical value", {
  for (tr in c("c", "ct")) {
    cv <- kpss_critical_values(tr)
    expect_identical(unname(vapply(cv, kpss_p_value, numeric(1), trend = tr)),
                     as.numeric(names(cv)))
  }
})
