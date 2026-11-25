# Shared argument checks
#
# Every user-supplied order, flag and significance level goes through these
# functions, so that a value is either accepted exactly as given or rejected
# with a clear message. Nothing is rounded, truncated or coerced silently.

# A non-negative whole number (lag order, maximum lag order, bandwidth,
# frequency, window size, number of replications). Returns it as an integer.
# Values above 1e6 are rejected: no lag order, frequency or window can be
# that large for a series that fits in memory, and the bound keeps later
# integer arithmetic (k + 3, 3s + 3, 10s) far from overflow (audit of 0.3.0).
check_order <- function(x, name) {
  ok <- (is.numeric(x) && !is.logical(x)) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x >= 0 && x == round(x)
  if (!ok) {
    stop("`", name, "` must be a single non-negative whole number ",
         "(e.g. 0, 1, 2); got ", describe_value(x), ".", call. = FALSE)
  }
  if (x > max_order) {
    stop("`", name, "` = ", format(x, big.mark = ",", scientific = FALSE),
         " is too large (at most ", format(max_order, big.mark = ",", scientific = FALSE), ").",
         call. = FALSE)
  }
  as.integer(x)
}

max_order <- 1e6

# A single TRUE or FALSE.
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE; got ", describe_value(x), ".",
         call. = FALSE)
  }
  x
}

# A significance level among those for which a decision can be made.
check_level <- function(level, available = c(0.01, 0.05, 0.10)) {
  ok <- is.numeric(level) && length(level) == 1L && !is.na(level) &&
    any(abs(available - level) < 1e-12)
  if (!ok) {
    stop("`level` must be one of ",
         paste(as.character(available), collapse = ", "),
         " (the levels with tabulated critical values); got ",
         describe_value(level), ".", call. = FALSE)
  }
  level
}

describe_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (length(x) != 1L) return(sprintf("a %s of length %d", class(x)[1], length(x)))
  if (is.character(x)) return(sprintf("\"%s\"", x))
  format(x)
}

# Information criterion for lag selection, in any case ("aic", "Maic"...):
# students type lower case. The default (the whole vector) gives "MAIC".
criterion_choices <- c("MAIC", "AIC", "BIC", "MBIC")
match_criterion <- function(criterion) {
  if (is.character(criterion) && length(criterion) == 1L && !is.na(criterion) &&
      nzchar(criterion)) {
    # full names or unambiguous prefixes, in any case (as match.arg allowed
    # before 0.3.2: "AI" is AIC; "M" is ambiguous)
    i <- pmatch(toupper(criterion), criterion_choices)
    if (!is.na(i)) return(criterion_choices[i])
  } else if (identical(criterion, criterion_choices)) {
    return("MAIC")
  }
  stop("`criterion` must be one of \"MAIC\", \"AIC\", \"BIC\", \"MBIC\" (any case, or an ",
       "unambiguous abbreviation); got ",
       describe_value(criterion), ".", call. = FALSE)
}

# A seed for set.seed(): a whole number within R's integer range (it takes no
# part in any arithmetic, so the 1e6 bound of check_order() does not apply).
check_seed <- function(seed) {
  ok <- is.numeric(seed) && !is.logical(seed) && length(seed) == 1L && !is.na(seed) &&
    is.finite(seed) && seed == round(seed) && abs(seed) <= .Machine$integer.max
  if (!ok) {
    stop("`seed` must be a single whole number between -", .Machine$integer.max,
         " and ", .Machine$integer.max, "; got ", describe_value(seed), ".", call. = FALSE)
  }
  as.integer(seed)
}

# Frequency of a series for the lag grids (residual checks, window sizes):
# an explicit `frequency` argument, else the series' own (validated like the
# argument, so a ts with an absurd frequency fails cleanly), else NA.
series_frequency <- function(frequency, series) {
  if (!is.null(frequency)) {
    s <- check_order(frequency, "frequency")
    if (s < 1L) stop("`frequency` must be at least 1.", call. = FALSE)
    return(s)
  }
  f <- if (is.null(series)) NA_real_ else series$frequency
  if (!is.finite(f) || f < 1) return(NA_integer_)
  if (f > max_order) {
    stop("The frequency of the series (", format(f, big.mark = ",", scientific = FALSE),
         ") is too large (at most ", format(max_order, big.mark = ",", scientific = FALSE),
         "); pass `frequency =` explicitly.", call. = FALSE)
  }
  as.integer(round(f))
}
