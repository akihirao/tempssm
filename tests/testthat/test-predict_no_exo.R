# tests/testthat/test-predict_no_exo

test_that(".predict_no_exo runs without error", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_train_mts <- .make_mts(y_train)

  res <- tempssm(y_train_mts)

  expect_no_error({
    .predict_no_exo(res$model, h = 12)
  })
})


test_that("prediction length matches h", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_train_mts <- .make_mts(y_train)

  res <- tempssm(y_train_mts)

  h <- 12
  pred <- .predict_no_exo(res$model, h)

  expect_equal(length(pred), h)
})


test_that("prediction output has valid type", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_train_mts <- .make_mts(y_train)

  res <- tempssm(y_train_mts)

  pred <- .predict_no_exo(res$model, h = 6)

  expect_true(
    is.numeric(pred) ||
      inherits(pred, "ts")
  )
})


test_that("different h produces different length", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_train_mts <- .make_mts(y_train)

  res <- tempssm(y_train_mts)

  pred1 <- .predict_no_exo(res$model, h = 3)
  pred2 <- .predict_no_exo(res$model, h = 6)

  expect_equal(length(pred1), 3)
  expect_equal(length(pred2), 6)
})


test_that("works with minimal horizon h = 1", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_train_mts <- .make_mts(y_train)

  res <- tempssm(y_train_mts)

  pred <- .predict_no_exo(res$model, h = 1)

  expect_equal(length(pred), 1)
})
