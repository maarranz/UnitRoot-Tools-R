# Deterministic specifications (spec: "Deterministic specifications")

trend_aliases <- c(
  n = "n", none = "n",
  c = "c", drift = "c", constant = "c", intercept = "c",
  ct = "ct", trend = "ct"
)

trend_labels <- c(
  n  = "none",
  c  = "constant",
  ct = "constant and linear trend"
)

#' Canonical deterministic specification
#'
#' Maps a deterministic specification, given by its canonical code or by
#' one of the aliases used in other packages, to the canonical code used
#' throughout `unitroottools`.
#'
#' | Code   | Meaning                    | Accepted aliases                                  |
#' |--------|----------------------------|---------------------------------------------------|
#' | `"n"`  | no constant, no trend      | `"none"`                                          |
#' | `"c"`  | constant                   | `"drift"` (urca), `"constant"`, `"intercept"` (bootUR) |
#' | `"ct"` | constant and linear trend  | `"trend"`                                         |
#'
#' @param trend A single character string.
#' @param allowed Canonical codes accepted in this context (e.g. KPSS
#'   accepts only `"c"` and `"ct"`).
#' @return The canonical code: `"n"`, `"c"` or `"ct"`.
#' @examples
#' match_trend("drift")
#' match_trend("ct")
#' @export
match_trend <- function(trend, allowed = c("n", "c", "ct")) {
  if (!is.character(trend) || length(trend) != 1L || is.na(trend)) {
    stop("`trend` must be a single character string.", call. = FALSE)
  }
  code <- unname(trend_aliases[tolower(trend)])
  if (is.na(code)) {
    stop(sprintf(
      "Unknown deterministic specification \"%s\". Use one of: %s.",
      trend, paste0('"', names(trend_aliases), '"', collapse = ", ")
    ), call. = FALSE)
  }
  if (!code %in% allowed) {
    stop(sprintf(
      "Deterministic specification \"%s\" is not available here. Allowed: %s.",
      code, paste0('"', allowed, '"', collapse = ", ")
    ), call. = FALSE)
  }
  code
}
