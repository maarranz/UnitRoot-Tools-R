# High-precision references shared with the Python implementation (spec,
# Decision D10): 18 ADF regressions (three series: a 120-observation series and
# two copies rescaled to order 1e-12, one sign-reversed; trends "c" and "ct";
# 0, 1 and 4 lags), computed by exact binary64-to-decimal elimination at 80
# and 120 digits, independently of any production code. The files are
# byte-identical copies of UnitRoot-Tools validation/scaled_svd/adf_inputs.json
# and validation/enders/scaled_references.json.

fixture_path <- function(file) {
  system.file("fixtures", "python", file, package = "unitroottools")
}

test_that("the reference files are the Python originals", {
  skip_if_not_installed("jsonlite")
  inputs <- fixture_path("adf_inputs.json")
  refs <- jsonlite::fromJSON(fixture_path("scaled_references.json"),
                             simplifyVector = FALSE)
  # MD5 of the copies (base R); their SHA-256 values are recorded in
  # inst/fixtures/python/README.md, and the input SHA-256 also in the
  # reference file itself.
  expect_identical(unname(tools::md5sum(inputs)), "3fbe15cd3e4535779df1c7a66fdca63f")
  expect_identical(unname(tools::md5sum(fixture_path("scaled_references.json"))),
                   "a37297ee740f20dfafae85dff2642d37")
  expect_identical(refs$input_sha256,
                   "0a5ef87c88547da466d37b1368e195037fe010de7a32a75e9234d8f2524ae8d4")
  expect_length(refs$cases, 18L)
})

test_that("ADF regressions, tau, phi and auxiliary t match the 120-digit references", {
  skip_if_not_installed("jsonlite")
  inputs <- jsonlite::fromJSON(fixture_path("adf_inputs.json"))$inputs
  refs <- jsonlite::fromJSON(fixture_path("scaled_references.json"),
                             simplifyVector = FALSE)$cases
  tol <- 1e-10                       # spec default (relative)
  rename <- function(nm) {
    nm <- sub("^Level\\.L1$", "y_lag1", nm)
    nm <- sub("^Diff\\.L", "dy_lag", nm)
    nm
  }
  for (case in refs) {
    y <- as.numeric(inputs[[case$key]])          # hexadecimal binary64 values
    k <- case$lags
    trend <- case$trend
    label <- sprintf("%s, trend %s, k = %d", case$key, trend, k)
    r <- adf_regression(y, trend, k, first = k + 2L)
    expect_identical(r$nobs, as.integer(case$nobs), label = label)

    u <- case$fits$unrestricted
    nm <- rename(unlist(u$names))
    est <- setNames(unlist(u$params), nm)
    tv <- setNames(unlist(u$tvalues), nm)
    # Python's trend restarts at 1 in the regression sample; R uses the
    # original index t. Only the constant (and its t ratio) depends on that.
    check <- if (trend == "ct") setdiff(nm, "const") else nm
    expect_equal(r$fit$coefficients[check], est[check], tolerance = tol, label = label)
    expect_equal(r$fit$t_values[check], tv[check], tolerance = tol, label = label)
    expect_equal(r$fit$rss, as.numeric(u$rss), tolerance = tol, label = label)
    expect_identical(as.integer(r$fit$df), as.integer(u$df_resid), label = label)
    expect_equal(unname(r$tau), tv[["y_lag1"]], tolerance = tol, label = label)

    phi_name <- if (trend == "ct") "phi3" else "phi1"
    expect_equal(r$phi[[phi_name]], case$phi, tolerance = tol, label = label)

    a <- case$fits$auxiliary
    st <- strategy_statistics(y, k, k + 2L, 0.05)
    aux <- if (trend == "ct") st$tests$aux_trend else st$tests$aux_const
    term <- if (trend == "ct") "trend" else "const"
    a_t <- setNames(unlist(a$tvalues), rename(unlist(a$names)))
    expect_equal(aux$statistic, a_t[[term]], tolerance = tol, label = label)
    expect_identical(as.integer(aux$df), as.integer(a$df_resid), label = label)
  }
})
