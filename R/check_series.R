# Input series (spec: "Input series")
#
# The input is a real, finite, one-dimensional, equally spaced series.
# Nothing is dropped, reordered or modified silently: anything ambiguous is
# an error with a message saying what to do.
#
# Accepted:
#   * numeric vectors, ts, zoo and xts (xts is a zoo subclass);
#   * one-column matrices and data frames;
#   * multi-column objects (mts, matrices, data frames, tibbles, tsibbles)
#     when `var` names the column to test, or when there is exactly one
#     numeric column (e.g. a tibble with a date column and a value column);
#   * tsibbles holding a single series (at most one key value), regular and
#     without implicit gaps.

check_series <- function(y, var = NULL, arg = "y") {
  if (!is.null(var) && (!is.character(var) || length(var) != 1L || is.na(var))) {
    stop("`var` must be a single column name (a character string).", call. = FALSE)
  }
  if (length(dim(y)) > 2L) {
    stop("`", arg, "` is an array with ", length(dim(y)), " dimensions; ",
         "a unit-root test needs a single series.", call. = FALSE)
  }

  if (inherits(y, "tbl_ts")) {
    y <- tsibble_series(y, var, arg)
  } else if (inherits(y, "zoo")) {
    if (!requireNamespace("zoo", quietly = TRUE)) {
      stop("Package 'zoo' is needed to use a zoo series as `", arg, "`.",
           call. = FALSE)
    }
    if (!zoo::is.regular(y, strict = TRUE)) {
      stop("`", arg, "` is an irregularly spaced zoo series; ",
           "unit-root tests need equally spaced observations.", call. = FALSE)
    }
    y <- pick_column(as.matrix(zoo::coredata(y)), var, arg)
  } else if (stats::is.ts(y) || is.matrix(y) || is.data.frame(y)) {
    if (is.data.frame(y)) check_date_order(y, arg)
    y <- pick_column(y, var, arg)
  } else if (!is.null(var)) {
    stop("`var` can only be used when `", arg, "` has columns.", call. = FALSE)
  }

  if (!is.null(dim(y)) && (length(dim(y)) > 2L || NCOL(y) != 1L)) {
    stop("The selected series has ", paste(dim(y), collapse = " x "),
         " values in a matrix or array; a unit-root test needs a single ",
         "column of observations.", call. = FALSE)
  }
  if (!is.numeric(y) || is.logical(y)) {
    stop("`", arg, "` must be numeric.", call. = FALSE)
  }
  bad <- which(!is.finite(y))
  if (length(bad) > 0L) {
    shown <- paste(utils::head(bad, 5L), collapse = ", ")
    stop("`", arg, "` has ", length(bad), " missing or non-finite value(s) ",
         "(position", if (length(bad) > 1L) "s", " ", shown,
         if (length(bad) > 5L) ", ...", "). Nothing is dropped automatically: ",
         "remove or replace them explicitly.", call. = FALSE)
  }
  as.numeric(unname(y))
}

# One column from a matrix-like object: the one named by `var`, the only
# column, or the only numeric column.
pick_column <- function(x, var, arg) {
  if (is.data.frame(x)) {
    cols <- names(x)
    numeric_cols <- cols[vapply(x, function(v) {
      is.numeric(v) && !is.logical(v) && is.null(dim(v))
    }, logical(1))]
    get <- function(nm) x[[nm]]
  } else {
    x <- as.matrix(x)
    cols <- colnames(x)
    numeric_cols <- if (is.numeric(x)) cols else character()
    get <- function(nm) x[, nm]
    if (NCOL(x) == 1L && is.null(var)) return(as.vector(x))
  }

  if (!is.null(var)) {
    if (is.null(cols) || !var %in% cols) {
      stop("`", arg, "` has no column named \"", var, "\".", call. = FALSE)
    }
    return(get(var))
  }
  if (length(cols) == 1L) return(get(cols))
  if (length(numeric_cols) == 1L) return(get(numeric_cols))
  stop("`", arg, "` has ", length(numeric_cols), " numeric columns",
       if (length(numeric_cols) > 0L)
         paste0(" (", paste0("\"", numeric_cols, "\"", collapse = ", "), ")"),
       "; say which one to test with `var`.", call. = FALSE)
}

# Data frames carry no frequency, so regular spacing cannot be checked. When
# there is exactly one date-time column, at least require it to be strictly
# increasing (no unsorted rows or duplicated dates).
check_date_order <- function(x, arg) {
  is_time <- vapply(x, function(v) inherits(v, c("Date", "POSIXt")), logical(1))
  if (sum(is_time) != 1L) return(invisible())
  d <- x[[which(is_time)]]
  if (anyNA(d) || any(diff(as.numeric(d)) <= 0)) {
    stop("The date column \"", names(x)[is_time], "\" of `", arg,
         "` is not strictly increasing (unsorted rows, repeated or missing ",
         "dates). Sort the data or fix the dates first.", call. = FALSE)
  }
  invisible()
}

tsibble_series <- function(x, var, arg) {
  if (!requireNamespace("tsibble", quietly = TRUE)) {
    stop("Package 'tsibble' is needed to use a tsibble as `", arg, "`.",
         call. = FALSE)
  }
  keys <- tsibble::key_vars(x)
  if (length(keys) > 0L && tsibble::n_keys(x) > 1L) {
    stop("`", arg, "` contains ", tsibble::n_keys(x), " series (key: ",
         paste(keys, collapse = ", "), "). Filter it to a single series first.",
         call. = FALSE)
  }
  if (!tsibble::is_regular(x)) {
    stop("`", arg, "` is an irregular tsibble; unit-root tests need equally ",
         "spaced observations.", call. = FALSE)
  }
  if (any(tsibble::has_gaps(x)$.gaps)) {
    stop("`", arg, "` has implicit gaps (missing time periods). Fill them ",
         "explicitly, e.g. with tsibble::fill_gaps(), and decide how to treat ",
         "the missing values.", call. = FALSE)
  }
  idx <- tsibble::index_var(x)
  x <- as.data.frame(x)
  x <- x[order(x[[idx]]), , drop = FALSE]
  measured <- setdiff(names(x), c(idx, keys))
  pick_column(x[, measured, drop = FALSE], var, arg)
}
