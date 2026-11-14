# Fixes found while writing the unit-root Labs (0.3.2)

test_that("subsetting a recursive_adf() result keeps it printable (item 1)", {
  set.seed(1)
  y <- stats::ts(cumsum(stats::rnorm(200)), start = c(1960, 1), frequency = 12)
  rr <- recursive_adf(y, lags = 1)
  rows <- rr[rr$p_value < 0.9, ]
  expect_s3_class(rows, "ur_recursive")
  expect_identical(attr(rows, "settings"), attr(rr, "settings"))
  expect_match(capture.output(print(rows)), "Recursive ADF test", all = FALSE)
  cols <- rr[, c("end", "statistic")]
  expect_identical(class(cols), "data.frame")
  expect_null(attr(cols, "settings"))
  expect_silent(capture.output(print(cols)))
  expect_type(rr$statistic, "double")                # column access unchanged
  expect_identical(rr[["lags"]], rr$lags)
  # objects that lost their settings some other way still print
  bare <- rr; attr(bare, "settings") <- NULL
  expect_silent(capture.output(print(bare)))
  expect_error(plot(bare), "no longer has the settings")
  empty <- rr[rr$statistic < -100, ]
  expect_silent(capture.output(print(empty)))
})

test_that("criterion is accepted in any case (item 2)", {
  set.seed(2)
  y <- cumsum(stats::rnorm(150))
  for (cr in c("aic", "Aic", "bic", "maic", "mbic")) {
    expect_identical(adf(y, criterion = cr)$lag_selection$method, toupper(cr))
  }
  expect_identical(select_lags(y, criterion = "bic")$criterion, "BIC")
  expect_identical(dfgls(y, criterion = "aic")$lag_selection$method, "AIC")
  expect_identical(ng_perron(y, criterion = "aic")$lag_selection$method, "AIC")
  expect_identical(ers_pt(y, criterion = "aic")$lag_selection$method, "AIC")
  expect_s3_class(ur_strategy(y, criterion = "aic"), "ur_strategy")
  expect_identical(attr(recursive_adf(y, criterion = "aic", reselect = TRUE), "settings")$lag_choice,
                   "AIC, selected in each window (kmax from the window size)")
  expect_identical(adf(y)$lag_selection$method, "MAIC")   # default unchanged
  expect_error(adf(y, criterion = "foo"), "any case")
  expect_error(adf(y, criterion = c("aic", "bic")), "must be one of")
  expect_error(adf(y, criterion = NA_character_), "must be one of")
})

test_that("printed samples show dates when the series has them (item 3)", {
  set.seed(3)
  v <- cumsum(stats::rnorm(200))
  m <- stats::ts(v, start = c(1960, 1), frequency = 12)
  line <- function(x, pat) grep(pat, capture.output(print(x)), value = TRUE)
  expect_match(line(adf(m, lags = 1), "Effective sample"), "N = 198, 1960-03 to 1976-08 \\(t = 3, ..., 200\\)")
  expect_match(line(ur_strategy(m, lags = 1), "Sample"), "1960-03 to 1976-08")
  expect_match(line(select_lags(m), "Common sample"), "to 1976-08 \\(t = ")
  expect_match(line(dfgls(m, lags = 1), "Effective sample"), "1960-03 to 1976-08")
  q <- stats::ts(v[1:80], start = c(1990, 1), frequency = 4)
  expect_match(line(adf(q, lags = 1), "Effective sample"), "1990 Q3 to 2009 Q4")
  a <- stats::ts(v[1:80], start = 1900)
  expect_match(line(adf(a, lags = 1), "Effective sample"), "1902 to 1979")
  d <- data.frame(date = seq(as.Date("2000-01-01"), by = "month", length.out = 120), x = v[1:120])
  expect_match(line(adf(d, lags = 1), "Effective sample"), "2000-03 to 2009-12")
  # a plain vector has no dates: indices only, as before
  expect_match(line(adf(v, lags = 1), "Effective sample"), "^  Effective sample: N = 198 \\(t = 3, ..., 200\\)$")
})
