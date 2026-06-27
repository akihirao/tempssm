# tests/testthat/test-ts_slice.R

test_that("ts_slice works correctly", {
  sliced <- .ts_slice(temp_ts_test, 3, 6)

  expect_length(sliced, 4)
  expect_identical(
    as.numeric(sliced),
    as.numeric(temp_ts_test[3:6])
  )
  expect_identical(stats::frequency(sliced), stats::frequency(temp_ts_test))
})


test_that("ts_slice preserves multivariate observations", {
  multi_ts <- ts(
    cbind(x = 1:8, y = 11:18),
    start = c(2000, 1),
    frequency = 4
  )

  sliced <- .ts_slice(multi_ts, 3, 6)

  expect_identical(NROW(sliced), 4L)
  expect_identical(NCOL(sliced), 2L)
  expect_identical(colnames(sliced), c("x", "y"))
  expect_identical(as.numeric(sliced[, "x"]), as.numeric(3:6))
  expect_identical(as.numeric(sliced[, "y"]), as.numeric(13:16))
})


test_that("ts_slice validates inputs and observation bounds", {
  multi_ts <- ts(
    cbind(x = 1:8, y = 11:18),
    frequency = 4
  )

  expect_null(.ts_slice(NULL, 1, 2))
  expect_error(
    .ts_slice(1:8, 1, 2),
    "Input must be a <ts> object"
  )
  expect_error(
    .ts_slice(multi_ts, 7, 10),
    "Invalid slice indices"
  )
  expect_error(
    .ts_slice(multi_ts, c(1, 2), 4),
    "i_start.*length one"
  )
  expect_error(
    .ts_slice(multi_ts, 1, "4"),
    "i_end.*numeric"
  )
  expect_error(
    .ts_slice(multi_ts, NA_real_, 4),
    "Invalid slice indices"
  )
  expect_error(
    .ts_slice(multi_ts, 1.5, 4),
    "Invalid slice indices"
  )
  expect_error(
    .ts_slice(multi_ts, 0, 4),
    "Invalid slice indices"
  )
  expect_error(
    .ts_slice(multi_ts, 5, 4),
    "Invalid slice indices"
  )
})
