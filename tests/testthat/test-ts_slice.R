# tests/testthat/test-ts_slice.R

test_that("ts_slice works correctly", {
  sliced <- tempssm:::.ts_slice(temp_ts_test, 3, 6)

  expect_length(sliced, 4)
})
