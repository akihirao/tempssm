# tests/testthat/test-check_na_ratio.R

test_that(".check_na_ratio warns when threshold exceeded", {
  x <- c(NA, NA, 1, 2) # 50% NA

  expect_warning(
    .check_na_ratio(x, threshold = 0.3, msg = "too many NA"),
    "too many NA"
  )
})


test_that(".check_na_ratio does not warn when below threshold", {
  x <- c(NA, 1, 2, 3) # 25% NA

  expect_no_warning(
    .check_na_ratio(x, threshold = 0.3, msg = "should not trigger")
  )
})


test_that("no warning when exactly equal to threshold", {
  x <- c(NA, NA, 1, 2) # 50%

  expect_no_warning(
    .check_na_ratio(x, threshold = 0.5, msg = "edge case")
  )
})


test_that("no warning when no missing values", {
  x <- c(1, 2, 3, 4)

  expect_no_warning(
    .check_na_ratio(x, threshold = 0.1, msg = "no NA")
  )
})


test_that("warning for all NA", {
  x <- c(NA, NA, NA)

  expect_warning(
    .check_na_ratio(x, threshold = 0.3, msg = "all NA"),
    "all NA"
  )
})
