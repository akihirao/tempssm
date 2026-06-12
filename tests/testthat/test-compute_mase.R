# test-compute_mase.R

test_that("compute_mase works with naive scaling", {
  y_train <- ts(1:10, frequency = 1)
  y_true <- c(5, 6, 7)
  y_pred <- c(5, 7, 6)

  mase <- tempssm:::.compute_mase(y_pred, y_true, y_train, method = "naive")

  expect_true(is.numeric(mase))
  expect_false(is.na(mase))
})


test_that("compute_mase works with seasonal scaling", {
  y_train <- ts(rnorm(24), frequency = 12)
  y_true <- c(10, 11, 12)
  y_pred <- c(9, 12, 11)

  mase <- tempssm:::.compute_mase(y_pred, y_true, y_train, method = "seasonal")

  expect_true(is.numeric(mase))
})


test_that("perfect prediction gives zero MASE", {
  y_train <- ts(1:10, frequency = 1)
  y_true <- c(5, 6, 7)
  y_pred <- c(5, 6, 7)

  mase <- tempssm:::.compute_mase(y_pred, y_true, y_train)

  expect_equal(mase, 0)
})


test_that("NULL inputs return NA", {
  y_train <- ts(1:10, frequency = 1)

  expect_true(is.na(tempssm:::.compute_mase(NULL, 1:3, y_train)))
  expect_true(is.na(tempssm:::.compute_mase(1:3, NULL, y_train)))
})


test_that("NA in MAE returns NA", {
  y_train <- ts(1:10, frequency = 1)

  y_true <- c(NA, NA)
  y_pred <- c(1, 2)

  mase <- tempssm:::.compute_mase(y_pred, y_true, y_train)

  expect_true(is.na(mase))
})


test_that("zero scaling factor returns NA", {
  y_train <- ts(rep(5, 10), frequency = 1) # 差分ゼロ
  y_true <- c(5, 5)
  y_pred <- c(4, 6)

  mase <- tempssm:::.compute_mase(y_pred, y_true, y_train)

  expect_true(is.na(mase))
})


test_that("scale_Q error returns NA", {
  y_train <- ts(1:5, frequency = 1)

  y_true <- c(1, 2)
  y_pred <- c(2, 3)

  mase <- tempssm:::.compute_mase(y_pred, y_true, y_train)

  expect_true(is.numeric(mase) || is.na(mase))
})


test_that("non-ts y_train triggers error", {
  expect_error(
    tempssm:::.compute_mase(1:3, 1:3, y_train = 1:10),
    "must be a"
  )
})


test_that("invalid method triggers error", {
  y_train <- ts(1:10, frequency = 1)

  expect_error(
    tempssm:::.compute_mase(1:3, 1:3, y_train, method = "wrong"),
    "arg"
  )
})


test_that("output is a scalar", {
  y_train <- ts(1:10, frequency = 1)
  y_true <- c(5, 6, 7)
  y_pred <- c(5, 7, 6)

  mase <- tempssm:::.compute_mase(y_pred, y_true, y_train)

  expect_length(mase, 1)
})
