test_that("MASE calculations match fixed hand-computed values", {
  y_train_naive <- ts(c(1, 2, 4, 7), frequency = 1)
  y_true <- c(10, 14)
  y_pred <- c(8, 17)

  naive_mase <- .compute_mase(
    y_pred,
    y_true,
    y_train_naive,
    method = "naive"
  )

  expect_equal(naive_mase, 1.25)

  y_train_seasonal <- ts(c(1, 10, 3, 16, 5, 22), frequency = 2)

  seasonal_mase <- .compute_mase(
    y_pred,
    y_true,
    y_train_seasonal,
    method = "seasonal"
  )

  expect_equal(seasonal_mase, 0.625)
})


test_that("temperature anomalies match fixed seasonal means", {
  temp_ts <- ts(
    c(10, 20, 12, 22),
    start = c(2000, 1),
    frequency = 2
  )

  anomaly <- compute_temp_anomaly(temp_ts)

  expect_identical(as.numeric(anomaly), c(-1, -1, 1, 1))
  expect_identical(stats::start(anomaly), c(2000, 1))
  expect_identical(stats::frequency(anomaly), 2)
})


test_that("monthly data-frame conversion matches a fixed regular ts", {
  temp_df <- data.frame(
    Date = as.Date(c("2001-01-01", "2001-02-01", "2001-03-01")),
    Temp = c(10.5, 11.5, 12.5)
  )

  temp_ts <- convert_monthly_df_to_ts(temp_df)

  expect_identical(as.numeric(temp_ts), c(10.5, 11.5, 12.5))
  expect_identical(stats::start(temp_ts), c(2001, 1))
  expect_identical(stats::frequency(temp_ts), 12)
})


test_that("rolling-origin split indices match fixed expected values", {
  temp_ts <- ts(seq_len(10), start = c(2000, 1), frequency = 4)

  folds <- ts_train_test_split(
    temp_ts,
    initial = 4,
    horizon = 2,
    step = 2
  )

  expect_length(folds, 3)
  expect_identical(as.numeric(folds[[1]]$train_ts), as.numeric(1:4))
  expect_identical(as.numeric(folds[[1]]$test_ts), as.numeric(5:6))
  expect_identical(as.numeric(folds[[2]]$train_ts), as.numeric(1:6))
  expect_identical(as.numeric(folds[[2]]$test_ts), as.numeric(7:8))
  expect_identical(as.numeric(folds[[3]]$train_ts), as.numeric(1:8))
  expect_identical(as.numeric(folds[[3]]$test_ts), as.numeric(9:10))
  expect_identical(folds[[1]]$train_idx, c(1, 4))
  expect_identical(folds[[1]]$test_idx, c(5, 6))
})
