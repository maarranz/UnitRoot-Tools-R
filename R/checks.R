# Shared argument checks
#
# Every user-supplied order, flag and significance level goes through these
# functions, so that a value is either accepted exactly as given or rejected
# with a clear message. Nothing is rounded, truncated or coerced silently.

# A non-negative whole number (lag order, maximum lag order, bandwidth).
# Returns it as an integer.
check_order <- function(x, name) {
  ok <- (is.numeric(x) && !is.logical(x)) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x >= 0 && x == round(x) &&
    x <= .Machine$integer.max
  if (!ok) {
    stop("`", name, "` must be a single non-negative whole number ",
         "(e.g. 0, 1, 2); got ", describe_value(x), ".", call. = FALSE)
  }
  as.integer(x)
}

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
