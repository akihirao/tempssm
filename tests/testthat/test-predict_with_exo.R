# tests/testthat/test-predict_with_exo


test_that(".predict_with_exo runs without error", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_test <- window(temp_ts_small, start = c(2003, 1))

  x_train <- window(exo_ts_small, end = c(2002, 12))
  x_test <- window(exo_ts_small, start = c(2003, 1))

  y_train_named <- set_ts_name(y_train, label = "Temp", quiet = TRUE)
  y_test_named <- set_ts_name(y_test, label = "Temp", quiet = TRUE)

  res <- tempssm(y_train_named, x_train)

  expect_no_error({
    .predict_with_exo(
      res,
      y_train_named,
      y_test_named,
      x_test,
      ar_order = 1,
      use_season = TRUE
    )
  })
})


# test_that("prediction length matches test length", {
#
#  y_train <- window(temp_ts_small, end = c(2002, 12))
#  y_test  <- window(temp_ts_small, start = c(2003, 1))
#
#  x_train <- window(exo_ts_small, end = c(2002, 12))
#  x_test  <- window(exo_ts_small, start = c(2003, 1))
#
#  y_train_named <- set_ts_name(y_train, label = "Temp", quiet = TRUE)
#  y_test_named <- set_ts_name(y_test, label = "Temp", quiet = TRUE)
#
#  res <- tempssm(y_train_named, x_train)
#
#  pred <- .predict_with_exo(
#    res,
#    y_train_named,
#    y_test_named,
#    x_test,
#    ar_order = 1,
#    use_season = FALSE
#  )
#
#  expect_identical(length(pred), length(y_test))
# })


test_that("prediction output has valid type", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_test <- window(temp_ts_small, start = c(2003, 1))

  x_train <- window(exo_ts_small, end = c(2002, 12))
  x_test <- window(exo_ts_small, start = c(2003, 1))

  y_train_named <- set_ts_name(y_train, label = "Temp", quiet = TRUE)
  y_test_named <- set_ts_name(y_test, label = "Temp", quiet = TRUE)

  res <- tempssm(y_train_named, x_train)

  pred <- .predict_with_exo(
    res,
    y_train_named,
    y_test_named,
    x_test,
    ar_order = 1,
    use_season = TRUE
  )

  expect_true(is.numeric(pred) || inherits(pred, "ts"))
})


test_that(".predict_with_exo reuses fitted tempssm object", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  y_test <- window(temp_ts_small, start = c(2003, 1))

  x_train <- window(exo_ts_small, end = c(2002, 12))
  x_test <- window(exo_ts_small, start = c(2003, 1))

  y_train_named <- set_ts_name(y_train, label = "Temp", quiet = TRUE)
  y_test_named <- set_ts_name(y_test, label = "Temp", quiet = TRUE)

  res <- tempssm(y_train_named, x_train)

  mockery::stub(
    .predict_with_exo,
    "tempssm",
    function(...) {
      stop("tempssm() should not be called", call. = FALSE)
    }
  )

  expect_no_error(
    .predict_with_exo(
      res,
      y_train_named,
      y_test_named,
      x_test,
      ar_order = 1,
      use_season = TRUE
    )
  )
})
