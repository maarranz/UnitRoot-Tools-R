# Specification strategies (spec: "Specification strategies")
#
# Structure:
#   strategy_statistics()  every statistic a strategy may need, all computed
#                          on one common sample with one lag order;
#   strategy_path()        the decision logic: a pure function of the
#                          rejection decisions (no data), so each branch can
#                          be tested directly and the whole procedure can be
#                          bootstrapped later;
#   ur_strategy()          the user-facing function.
#
# Notation (spec): b is the coefficient on y_{t-1}, gamma the trend
# coefficient, a the constant. Enders writes gamma for b and a2 for gamma.

strategy_names <- c(
  original = "Original (Dolado-Jenkinson-Sosvilla-Rivero / Enders)",
  conservative = "Conservative",
  hamilton = "Hamilton (rule-based)"
)

# Per-test descriptions used in the history.
strategy_tests <- list(
  tau3      = list(stage = "T", regression = "ct",     null = "b = 0: unit root (ct regression)"),
  phi3      = list(stage = "T", regression = "ct",     null = "b = gamma = 0: unit root and no trend"),
  aux_trend = list(stage = "T", regression = "aux_ct", null = "gamma = 0: no trend, given b = 0"),
  cond_ct   = list(stage = "T", regression = "ct",     null = "b = 0: unit root, trend retained"),
  tau2      = list(stage = "C", regression = "c",      null = "b = 0: unit root (c regression)"),
  phi1      = list(stage = "C", regression = "c",      null = "a = b = 0: unit root and no constant"),
  aux_const = list(stage = "C", regression = "aux_c",  null = "a = 0: no drift, given b = 0"),
  cond_c    = list(stage = "C", regression = "c",      null = "b = 0: unit root, constant retained"),
  tau1      = list(stage = "N", regression = "n",      null = "b = 0: unit root (n regression)")
)

rejection_types <- c(
  tau3 = "unit_root", tau2 = "unit_root", tau1 = "unit_root",
  cond_ct = "unit_root", cond_c = "unit_root",
  phi3 = "joint_restriction", phi1 = "joint_restriction",
  aux_trend = "deterministic", aux_const = "deterministic"
)

# ---- statistics ----------------------------------------------------------------

strategy_statistics <- function(y, k, first, level) {
  lv <- format_level(level)
  r <- list(ct = adf_regression(y, "ct", k, first),
            c  = adf_regression(y, "c",  k, first),
            n  = adf_regression(y, "n",  k, first))
  nobs <- r$ct$nobs

  tau_entry <- function(trend) {
    cv <- tau_critical_values(trend, nobs)[[lv]]
    stat <- unname(r[[trend]]$tau)
    list(statistic = stat, critical_value = cv, reject = stat < cv,
         p_value = tau_p_value(stat, trend, nobs),
         calibration = sprintf("%s, N = %d", mackinnon_source[[trend]], nobs))
  }
  phi_entry <- function(trend, name) {
    stat <- r[[trend]]$phi[[name]]
    T_eq <- nobs + k + 1L
    p <- phi_calibration(name, stat, T_eq, k)
    cv <- p[[lv]]
    list(statistic = stat, critical_value = cv, reject = stat > cv,
         p_value = attr(p, "p_value"),
         calibration = sprintf("%s, T = %d, k = %d", phi_source, T_eq, k))
  }
  # t test of a deterministic coefficient with the unit root imposed (b = 0):
  # drop y_lag1 from the design; no I(1) regressors remain, so Student t.
  aux_entry <- function(trend, coef_name) {
    d <- r[[trend]]$design
    keep <- setdiff(colnames(d$X), "y_lag1")
    f <- qr_fit(d$X[, keep, drop = FALSE], d$response)
    stat <- f$t_values[[coef_name]]
    cv <- stats::qt(1 - level / 2, f$df)
    list(statistic = stat, critical_value = cv, reject = abs(stat) > cv,
         p_value = 2 * stats::pt(-abs(stat), f$df),
         calibration = sprintf("Student t, df = %d, two-sided", f$df), df = f$df)
  }
  # conditional t test of b in the regression with the retained deterministic
  # term (Enders' Rule 2): same t ratio as tau, Student t reference, lower tail.
  cond_entry <- function(trend) {
    f <- r[[trend]]$fit
    stat <- f$t_values[["y_lag1"]]
    cv <- stats::qt(level, f$df)
    list(statistic = stat, critical_value = cv, reject = stat < cv,
         p_value = stats::pt(stat, f$df),
         calibration = sprintf("Student t, df = %d, lower tail", f$df), df = f$df)
  }

  list(
    regressions = r,
    tests = list(
      tau3 = tau_entry("ct"), phi3 = phi_entry("ct", "phi3"),
      aux_trend = aux_entry("ct", "trend"), cond_ct = cond_entry("ct"),
      tau2 = tau_entry("c"), phi1 = phi_entry("c", "phi1"),
      aux_const = aux_entry("c", "const"), cond_c = cond_entry("c"),
      tau1 = tau_entry("n")
    )
  )
}

format_level <- function(level) {
  lv <- c("0.01", "0.05", "0.1")[match(TRUE, abs(c(0.01, 0.05, 0.10) - level) < 1e-12)]
  if (is.na(lv)) {
    stop("Specification strategies use tabulated critical values: `level` must be ",
         "0.01, 0.05 or 0.10.", call. = FALSE)
  }
  lv
}

# ---- decision logic (pure) -------------------------------------------------------
#
# `reject` is a named logical vector with one element per test in
# strategy_tests. Returns the visited steps with the next action and reason,
# and the terminal test, specification, decision and conclusion.

strategy_path <- function(reject, strategy, allow_no_constant = FALSE,
                          trending = NULL) {
  steps <- list()
  add <- function(test, next_action, reason) {
    steps[[length(steps) + 1L]] <<- data.frame(
      test = test, next_action = next_action, reason = reason,
      stringsAsFactors = FALSE)
  }
  finish <- function(terminal, specification, conclusion) {
    list(steps = do.call(rbind, steps), terminal = terminal,
         specification = specification,
         unit_root_rejected = unname(reject[[terminal]]),
         conclusion = conclusion)
  }

  if (strategy == "hamilton") {
    if (isTRUE(trending)) {
      add("tau3", "stop", "specification declared by user: trending series")
      return(finish("tau3", "ct", if (reject[["tau3"]])
        "no unit root: stationary around a linear trend" else
        "unit root (with drift)"))
    }
    add("tau2", "stop", "specification declared by user: non-trending series")
    return(finish("tau2", "c", if (reject[["tau2"]])
      "no unit root: stationary around a constant level" else
      "unit root (no drift)"))
  }

  # Stage T: constant and trend
  if (reject[["tau3"]]) {
    add("tau3", "stop", "unit root rejected")
    return(finish("tau3", "ct", "no unit root: stationary around a linear trend"))
  }
  add("tau3", "test phi3", "unit root not rejected; check whether the trend belongs")
  stage_c <- TRUE
  if (!reject[["phi3"]]) {
    add("phi3", "go to stage C", "b = gamma = 0 not rejected: drop the trend")
  } else {
    add("phi3", "test trend, b = 0 imposed",
        "b = gamma = 0 rejected: the trend may belong")
    if (!reject[["aux_trend"]]) {
      add("aux_trend", "go to stage C", "trend not significant: drop it")
    } else {
      add("aux_trend", "conditional t test of b",
          "trend significant: keep it; t distribution applies to b (Enders, Rule 2)")
      add("cond_ct", "stop", "terminal test with the trend retained")
      return(finish("cond_ct", "ct", if (reject[["cond_ct"]])
        "no unit root: stationary around a linear trend" else
        "unit root with drift and a trend in the differences (quadratic trend in levels)"))
    }
  }

  # Stage C: constant only
  if (reject[["tau2"]]) {
    add("tau2", "stop", "unit root rejected")
    return(finish("tau2", "c", "no unit root: stationary around a non-zero mean"))
  }
  if (strategy == "conservative" && !allow_no_constant) {
    add("tau2", "stop",
        "unit root not rejected; constant retained (no-constant model not considered)")
    return(finish("tau2", "c", "unit root (constant retained)"))
  }
  add("tau2", "test phi1", "unit root not rejected; check whether the constant belongs")
  if (!reject[["phi1"]]) {
    add("phi1", "go to stage N", "a = b = 0 not rejected: drop the constant")
  } else {
    add("phi1", "test constant, b = 0 imposed",
        "a = b = 0 rejected: the constant may belong")
    if (!reject[["aux_const"]]) {
      add("aux_const", "go to stage N", "drift not significant: drop the constant")
    } else {
      add("aux_const", "conditional t test of b",
          "drift significant: keep the constant; t distribution applies to b (Enders, Rule 2)")
      add("cond_c", "stop", "terminal test with the constant retained")
      return(finish("cond_c", "c", if (reject[["cond_c"]])
        "no unit root: stationary around a non-zero mean" else
        "unit root with drift"))
    }
  }

  # Stage N: no deterministic terms
  add("tau1", "stop", "terminal test without deterministic terms")
  finish("tau1", "n", if (reject[["tau1"]])
    "no unit root: stationary with zero mean" else "unit root without drift")
}

# ---- user function ------------------------------------------------------------

#' Specification strategies for unit-root testing
#'
#' Chooses the deterministic terms of the ADF regression and tests for a unit
#' root, recording every step of the decision.
#'
#' * `"original"`: the Dolado-Jenkinson-Sosvilla-Rivero sequence as set out
#'   in Enders (Supplementary Manual, Section 4.4): \eqn{\tau_3}, \eqn{\phi_3},
#'   a t test of the trend with the unit root imposed, and a conditional t test
#'   of the lagged level; then the same with a constant only (\eqn{\tau_2},
#'   \eqn{\phi_1}, ...); finally \eqn{\tau_1} without deterministic terms.
#' * `"conservative"`: as the original, but when \eqn{\tau_2} does not reject
#'   the constant is retained and the procedure stops. With
#'   `allow_no_constant = TRUE` the original constant-stage steps are
#'   evaluated. This is an experimental proposal: its overall size has not
#'   been established.
#' * `"hamilton"`: Hamilton's (1994, ch. 17) rule. The user declares whether
#'   the series is trending (`trending = TRUE`: \eqn{\tau_3} in the regression
#'   with constant and trend) or not (`trending = FALSE`: \eqn{\tau_2}). There
#'   is no default and no data-driven choice.
#'
#' All regressions use the same lag order and the same sample. The lag order
#' is selected once, under constant and trend (under the declared
#' specification for Hamilton's rule), unless fixed with `lags`; the
#' regressions then use the maximal sample for that order (or the common
#' selection sample with `final_sample = "common"`).
#'
#' The history distinguishes rejections of the unit-root null from
#' rejections that only decide the path of the procedure (joint \eqn{\phi}
#' restrictions and deterministic coefficients). No p-value is reported for
#' the procedure as a whole: its overall size is not known.
#'
#' @inheritParams adf
#' @param strategy `"original"`, `"conservative"` or `"hamilton"`.
#' @param trending For `strategy = "hamilton"`: `TRUE` if the series is
#'   trending, `FALSE` if not. Required; ignored by the other strategies.
#' @param allow_no_constant For `strategy = "conservative"`: allow the
#'   procedure to evaluate the model without a constant.
#' @param level Significance level for every test in the sequence: 0.01,
#'   0.05 or 0.10.
#'
#' @return An object of class `ur_strategy`: a list with the `strategy`,
#'   `settings`, chosen `specification`, `terminal` test (a `ur_test`),
#'   `unit_root_rejected`, `conclusion`, `lags`, `nobs`, `sample`,
#'   `lag_selection` and `history` (a data frame, one row per step).
#' @references
#' Dolado, J. J., Jenkinson, T. and Sosvilla-Rivero, S. (1990). Cointegration
#' and unit roots. *Journal of Economic Surveys*, 4, 249-273.
#'
#' Enders, W. *Supplementary Manual to accompany Applied Econometric Time
#' Series*, 4th ed., Section 4.4.
#'
#' Hamilton, J. D. (1994). *Time Series Analysis*. Princeton University Press,
#' ch. 17.
#' @examples
#' set.seed(1)
#' y <- cumsum(0.2 + rnorm(200))
#' ur_strategy(y, "original")
#' ur_strategy(y, "conservative")
#' ur_strategy(y, "hamilton", trending = TRUE)
#' @export
ur_strategy <- function(y, strategy = c("original", "conservative", "hamilton"),
                        trending = NULL, allow_no_constant = FALSE,
                        lags = NULL, max_lags = NULL,
                        criterion = c("MAIC", "AIC", "BIC", "MBIC"),
                        final_sample = c("maximal", "common"),
                        level = 0.05, var = NULL) {
  call <- match.call()
  strategy <- match.arg(strategy)
  criterion <- match.arg(criterion)
  final_sample <- match.arg(final_sample)
  check_level(level)
  check_flag(allow_no_constant, "allow_no_constant")
  y <- check_series(y, var = var)

  if (strategy == "hamilton") {
    if (is.null(trending) || !is.logical(trending) || length(trending) != 1L ||
        is.na(trending)) {
      stop("Hamilton's rule needs your judgement: set `trending = TRUE` for a ",
           "trending series (e.g. log GDP) or `trending = FALSE` for a ",
           "non-trending one (e.g. an interest rate). There is no default.",
           call. = FALSE)
    }
    selection_trend <- if (trending) "ct" else "c"
  } else {
    if (!is.null(trending)) {
      warning("`trending` is only used by strategy = \"hamilton\"; ignored.",
              call. = FALSE)
    }
    selection_trend <- "ct"
  }
  if (!isFALSE(allow_no_constant) && strategy != "conservative") {
    warning("`allow_no_constant` is only used by strategy = \"conservative\"; ignored.",
            call. = FALSE)
  }

  if (!is.null(lags)) {
    if (!is.null(max_lags)) {
      stop("Give either `lags` (fixed order) or `max_lags` (selection), not both.",
           call. = FALSE)
    }
    k <- check_order(lags, "lags")
    first <- k + 2L
    selection <- NULL
  } else {
    sel <- lag_table(y, selection_trend, max_lags, criterion)
    k <- sel$lags
    first <- if (final_sample == "common") sel$kmax + 2L else k + 2L
    selection <- list(method = criterion, kmax = sel$kmax,
                      detrending = sel$detrending, table = sel$table,
                      final_sample = final_sample, selected_under = selection_trend,
                      notes = sel$notes)
  }

  nobs <- length(y) - first + 1L
  if (nobs < 25L) {
    stop("Specification strategies need at least 25 effective observations ",
         "(here N = ", nobs, ").", call. = FALSE)
  }

  st <- strategy_statistics(y, k, first, level)
  reject <- vapply(st$tests, `[[`, logical(1), "reject")
  path <- strategy_path(reject, strategy, allow_no_constant, trending)

  history <- strategy_history(path$steps, st$tests, level)
  terminal <- terminal_test(path$terminal, st, selection, level, call)

  structure(list(
    strategy = strategy,
    settings = list(trending = trending, allow_no_constant = allow_no_constant,
                    level = level, criterion = if (is.null(lags)) criterion else NA),
    specification = path$specification,
    terminal = terminal,
    unit_root_rejected = path$unit_root_rejected,
    conclusion = path$conclusion,
    lags = k, nobs = nobs, sample = c(first, length(y)),
    lag_selection = selection,
    history = history,
    engine = engine_internal(),
    call = call
  ), class = "ur_strategy")
}

strategy_history <- function(steps, tests, level) {
  rows <- lapply(seq_len(nrow(steps)), function(i) {
    nm <- steps$test[i]
    t <- tests[[nm]]
    info <- strategy_tests[[nm]]
    data.frame(
      step = i, stage = info$stage, regression = info$regression, test = nm,
      null = info$null, statistic = t$statistic,
      critical_value = t$critical_value, p_value = t$p_value, level = level,
      calibration = t$calibration, reject = t$reject,
      rejection_type = if (t$reject) rejection_types[[nm]] else "none",
      next_action = steps$next_action[i], reason = steps$reason[i],
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

terminal_test <- function(name, st, selection, level, call) {
  tau_trend <- c(tau3 = "ct", tau2 = "c", tau1 = "n")
  if (name %in% names(tau_trend)) {
    return(adf_result(st$regressions[[tau_trend[[name]]]], selection, level, call))
  }
  trend <- if (name == "cond_ct") "ct" else "c"
  r <- st$regressions[[trend]]
  f <- r$fit
  stat <- f$t_values[["y_lag1"]]
  cv <- stats::qt(c(0.01, 0.05, 0.10), f$df)
  names(cv) <- c("0.01", "0.05", "0.1")
  coef_table <- cbind(Estimate = f$coefficients, `Std. Error` = f$se,
                      `t value` = f$t_values)
  new_ur_test(
    method = "Conditional t test of the lagged level",
    test = "ADF", statistic_name = "t_b",
    statistic = stat, null = "unit root", alternative = adf_alternatives[[trend]],
    trend = trend, lags = r$k, lag_selection = selection,
    nobs = r$nobs, sample = unname(r$sample),
    critical_values = cv, tail = "lower", p_value = stats::pt(stat, f$df),
    calibration = list(
      method = "conventional_t",
      source = sprintf("Student t with %d df", f$df),
      note = paste("valid only because the deterministic term was retained",
                   "(Enders, Rule 2); conditional on the earlier steps")),
    level = level, engine = engine_internal(),
    regression = list(coefficients = coef_table, rss = f$rss, df = f$df,
                      sigma2 = f$sigma2),
    call = call)
}

# ---- methods ------------------------------------------------------------------

test_labels <- c(tau3 = "tau3", phi3 = "phi3", aux_trend = "t(trend)",
                 cond_ct = "t(b) cond.", tau2 = "tau2", phi1 = "phi1",
                 aux_const = "t(const)", cond_c = "t(b) cond.", tau1 = "tau1")

#' @export
print.ur_strategy <- function(x, ...) {
  cat(sprintf("Specification strategy: %s\n\n", strategy_names[[x$strategy]]))
  ls <- x$lag_selection
  if (is.null(ls)) {
    cat(sprintf("  Lags: %d (fixed)\n", x$lags))
  } else {
    cat(sprintf("  Lags: %d (%s, kmax = %d, selected under \"%s\")\n",
                x$lags, ls$method, as.integer(ls$kmax), ls$selected_under))
    for (note in ls$notes) cat("  Note: ", note, "\n", sep = "")
  }
  cat(sprintf("  Sample (all regressions): N = %d (t = %d, ..., %d)   Level: %s\n\n",
              x$nobs, x$sample[1], x$sample[2], fmt_level(x$settings$level)))
  h <- x$history
  cat(sprintf("  %-4s %-6s %-10s %9s %9s %8s  %-10s  %s\n", "Step", "Regr.",
              "Test", "Statistic", "Crit.val.", "p-value", "Decision", "Next"))
  for (i in seq_len(nrow(h))) {
    pv <- if (is.na(h$p_value[i])) "-" else fmt_pvalue(h$p_value[i])
    cat(sprintf("  %-4d %-6s %-10s %9s %9s %8s  %-10s  %s\n", h$step[i],
                h$regression[i], test_labels[[h$test[i]]],
                fmt_num(h$statistic[i]), fmt_num(h$critical_value[i]), pv,
                if (h$reject[i]) "reject" else "not reject", h$next_action[i]))
  }
  cat(sprintf("\n  Specification: %s\n", trend_labels[[x$specification]]))
  cat(sprintf("  Conclusion:    %s\n", x$conclusion))
  cat("\n  Note: p-values are for each step on its own (none for the phi tables);\n",
      "  no p-value or size is reported for the sequence as a whole.\n", sep = "")
  invisible(x)
}

#' @export
summary.ur_strategy <- function(object, ...) {
  structure(list(result = object), class = "summary.ur_strategy")
}

#' @export
print.summary.ur_strategy <- function(x, ...) {
  r <- x$result
  print(r)
  cat("\nSteps\n")
  h <- r$history
  for (i in seq_len(nrow(h))) {
    cat(sprintf("  %d. %s   H0: %s\n", h$step[i], test_labels[[h$test[i]]], h$null[i]))
    cat(sprintf("     statistic %s, critical value %s%s [%s]: %s",
                fmt_num(h$statistic[i]), fmt_num(h$critical_value[i]),
                if (is.na(h$p_value[i])) "" else paste0(", p-value ", fmt_pvalue(h$p_value[i])),
                h$calibration[i], if (h$reject[i]) "reject" else "do not reject"))
    if (h$reject[i]) cat(sprintf(" (%s)", gsub("_", " ", h$rejection_type[i])))
    cat(sprintf("\n     -> %s: %s\n", h$next_action[i], h$reason[i]))
  }
  cat("\nTerminal test\n\n")
  print(r$terminal)
  invisible(x)
}

#' @export
as.data.frame.ur_strategy <- function(x, ...) x$history
