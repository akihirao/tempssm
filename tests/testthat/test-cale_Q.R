# test-cale_Q.R

test_that("scale_Q works for naive method", {
  ts_data <- ts(c(1, 2, 3, 5))

  Q <- .scale_Q(ts_data, method = "naive")

  expected <- mean(abs(diff(c(1, 2, 3, 5))))

  expect_identical(Q, expected)
})


test_that("scale_Q works for seasonal method", {
  ts_data <- ts(c(1, 2, 3, 4, 2, 3, 4, 5),
    frequency = 4
  )

  Q <- .scale_Q(ts_data, method = "seasonal")

  y <- as.numeric(ts_data)
  m <- 4
  n <- length(y)

  expected <- mean(abs(y[(m + 1):n] - y[1:(n - m)]))

  expect_identical(Q, expected)
})


test_that("scale_Q handles NA values with na.rm", {
  ts_data <- ts(c(1, NA, 3, 5))

  Q <- .scale_Q(ts_data, method = "naive")

  expect_type(Q, "double")
})


test_that("non-ts input triggers error", {
  expect_error(
    .scale_Q(1:10),
    "must be a"
  )
})


test_that("naive method requires at least 2 observations", {
  ts_data <- ts(1)

  expect_error(
    .scale_Q(ts_data, method = "naive"),
    "At least two observations"
  )
})


test_that("seasonal method requires frequency > 1", {
  ts_data <- ts(1:10, frequency = 1)

  expect_error(
    .scale_Q(ts_data, method = "seasonal"),
    "frequency > 1"
  )
})


test_that("seasonal method requires enough observations", {
  ts_data <- ts(1:4, frequency = 4)

  expect_error(
    .scale_Q(ts_data, method = "seasonal"),
    "too short"
  )
})


test_that("invalid method triggers error", {
  ts_data <- ts(1:10)

  expect_error(
    .scale_Q(ts_data, method = "invalid"),
    "arg"
  )
})


test_that("constant series returns zero Q", {
  ts_data <- ts(rep(5, 10), frequency = 1)

  Q <- .scale_Q(ts_data, method = "naive")

  expect_identical(Q, 0)
})


test_that("scale_Q returns scalar", {
  ts_data <- ts(1:10)

  Q <- .scale_Q(ts_data)

  expect_length(Q, 1)
})
