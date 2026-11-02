test_that("canonical codes are returned unchanged", {
  expect_identical(match_trend("n"), "n")
  expect_identical(match_trend("c"), "c")
  expect_identical(match_trend("ct"), "ct")
})

test_that("aliases from urca and bootUR are mapped", {
  expect_identical(match_trend("none"), "n")
  expect_identical(match_trend("drift"), "c")
  expect_identical(match_trend("constant"), "c")
  expect_identical(match_trend("intercept"), "c")
  expect_identical(match_trend("trend"), "ct")
  expect_identical(match_trend("Trend"), "ct")
})

test_that("unknown or disallowed specifications are errors", {
  expect_error(match_trend("ctt"), "Unknown deterministic specification")
  expect_error(match_trend("n", allowed = c("c", "ct")), "not available here")
  expect_error(match_trend(c("c", "ct")), "single character string")
  expect_error(match_trend(NA_character_), "single character string")
  expect_error(match_trend(1), "single character string")
})
