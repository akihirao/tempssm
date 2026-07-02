# tests/testthat/test-predict_with_exo.R

test_that(".predict_with_exo returns forecasts for future covariates", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  x_train <- window(exo_ts_small, end = c(2002, 12))
  x_test <- window(exo_ts_small, start = c(2003, 1))
  y_train <- set_ts_name(y_train, label = "Temp", quiet = TRUE)
  res <- tempssm(y_train, x_train)

  pred <- .predict_with_exo(res, x_test)

  expect_s3_class(pred, "ts")
  expect_length(pred, NROW(x_test))
})


test_that(".predict_with_exo supports prediction intervals", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  x_train <- window(exo_ts_small, end = c(2002, 12))
  x_test <- window(exo_ts_small, start = c(2003, 1))
  y_train <- set_ts_name(y_train, label = "Temp", quiet = TRUE)
  res <- tempssm(y_train, x_train)

  pred <- .predict_with_exo(
    res,
    x_test,
    interval = "prediction",
    level = 0.9
  )

  expect_s3_class(pred, "mts")
  expect_identical(colnames(pred), c("fit", "lwr", "upr"))
})


test_that(".predict_with_exo reuses the fitted tempssm object", {
  y_train <- window(temp_ts_small, end = c(2002, 12))
  x_train <- window(exo_ts_small, end = c(2002, 12))
  x_test <- window(exo_ts_small, start = c(2003, 1))
  y_train <- set_ts_name(y_train, label = "Temp", quiet = TRUE)
  res <- tempssm(y_train, x_train)

  mockery::stub(
    .predict_with_exo,
    "tempssm",
    function(...) {
      stop("tempssm() should not be called", call. = FALSE)
    }
  )

  expect_no_error(.predict_with_exo(res, x_test))
})
