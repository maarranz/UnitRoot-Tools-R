test_that("numeric vectors and ts objects are accepted and returned as plain numeric", {
  expect_identical(check_series(c(1, 2.5, 3)), c(1, 2.5, 3))
  expect_identical(check_series(1:3), c(1, 2, 3))
  y <- ts(c(1, 4, 2, 8), start = c(2000, 1), frequency = 4)
  expect_identical(check_series(y), c(1, 4, 2, 8))
  expect_identical(check_series(matrix(1:4, ncol = 1)), c(1, 2, 3, 4))
})

test_that("missing and non-finite values are errors, with positions", {
  expect_error(check_series(c(1, NA, 3)), "1 missing or non-finite value\\(s\\) \\(position 2\\)")
  expect_error(check_series(c(1, Inf, NaN, 4)), "positions 2, 3")
})

test_that("non-numeric inputs are errors", {
  expect_error(check_series(c("a", "b")), "must be numeric")
  expect_error(check_series(c(TRUE, FALSE)), "must be numeric")
  expect_error(check_series(1:3, var = 2), "single column name")
  expect_error(check_series(1:3, var = "x"), "only be used when")
})

test_that("multivariate ts and matrices need `var`", {
  m <- ts(cbind(gdp = 1:5, cpi = 11:15), start = 2000)
  expect_error(check_series(m), "2 numeric columns \\(\"gdp\", \"cpi\"\\); say which one")
  expect_identical(check_series(m, var = "cpi"), as.numeric(11:15))
  expect_error(check_series(m, var = "unemp"), "no column named \"unemp\"")
})

test_that("data frames: the only numeric column is used, or the one in `var`", {
  df <- data.frame(date = as.Date("2020-01-01") + 0:3, value = c(1, 3, 2, 5))
  expect_identical(check_series(df), c(1, 3, 2, 5))
  df2 <- data.frame(date = as.Date("2020-01-01") + 0:3, a = 1:4, b = c(2, 4, 6, 9))
  expect_error(check_series(df2), "say which one to test with `var`")
  expect_identical(check_series(df2, var = "b"), c(2, 4, 6, 9))
})

test_that("data frames with an unsorted or repeated date column are errors", {
  df <- data.frame(date = as.Date("2020-01-01") + c(0, 2, 1), value = 1:3)
  expect_error(check_series(df), "not strictly increasing")
  df <- data.frame(date = as.Date("2020-01-01") + c(0, 1, 1), value = 1:3)
  expect_error(check_series(df), "not strictly increasing")
})

test_that("tibbles are handled like data frames", {
  skip_if_not_installed("tibble")
  tb <- tibble::tibble(date = as.Date("2020-01-01") + 0:3, y = c(5, 6, 4, 7))
  expect_identical(check_series(tb), c(5, 6, 4, 7))
  tb2 <- tibble::tibble(date = as.Date("2020-01-01") + 0:3, a = 1:4, b = 4:1)
  expect_identical(check_series(tb2, var = "b"), c(4, 3, 2, 1))
})

test_that("zoo and xts series must be regular", {
  skip_if_not_installed("zoo")
  z <- zoo::zoo(c(1, 2, 3, 4), as.Date("2020-01-01") + c(0, 1, 2, 3))
  expect_identical(check_series(z), c(1, 2, 3, 4))
  zi <- zoo::zoo(c(1, 2, 3, 4), as.Date("2020-01-01") + c(0, 1, 5, 6))
  expect_error(check_series(zi), "irregularly spaced")
  zm <- zoo::zoo(cbind(a = 1:4, b = 5:8), as.Date("2020-01-01") + 0:3)
  expect_identical(check_series(zm, var = "b"), c(5, 6, 7, 8))
  skip_if_not_installed("xts")
  x <- xts::xts(c(2, 4, 3, 5), as.Date("2020-01-01") + 0:3)
  expect_identical(check_series(x), c(2, 4, 3, 5))
})

test_that("tsibbles: single regular series without gaps", {
  skip_if_not_installed("tsibble")
  tsb <- tsibble::tsibble(
    month = tsibble::yearmonth("2020-01") + 0:5,
    gdp = c(1, 2, 4, 3, 5, 6),
    index = month
  )
  expect_identical(check_series(tsb), c(1, 2, 4, 3, 5, 6))

  gaps <- tsb[-3, ]
  expect_error(check_series(gaps), "implicit gaps")

  two <- tsibble::tsibble(
    month = rep(tsibble::yearmonth("2020-01") + 0:2, 2),
    country = rep(c("ES", "FR"), each = 3),
    gdp = 1:6,
    index = month, key = country
  )
  expect_error(check_series(two), "contains 2 series \\(key: country\\)")

  wide <- tsibble::tsibble(
    month = tsibble::yearmonth("2020-01") + 0:3,
    gdp = 1:4, cpi = c(2, 3, 5, 4),
    index = month
  )
  expect_error(check_series(wide), "say which one to test with `var`")
  expect_identical(check_series(wide, var = "cpi"), c(2, 3, 5, 4))
})
