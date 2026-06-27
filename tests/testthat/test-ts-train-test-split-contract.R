test_that("rolling bounds support step values different from horizon", {
  temp_data <- ts(seq_len(14), start = c(2000, 1), frequency = 4)

  folds <- ts_train_test_split(
    temp_data,
    initial = 5,
    horizon = 3,
    step = 2
  )
  bounds <- vapply(
    folds,
    function(fold) c(fold$train_idx, fold$test_idx),
    numeric(4)
  )
  expected <- cbind(
    c(1, 5, 6, 8),
    c(1, 7, 8, 10),
    c(1, 9, 10, 12),
    c(1, 11, 12, 14)
  )

  expect_identical(bounds, expected)
})


test_that("fixed windows preserve their exact length and bounds", {
  temp_data <- ts(seq_len(14), start = c(2000, 1), frequency = 4)

  folds <- ts_train_test_split(
    temp_data,
    initial = 5,
    horizon = 3,
    step = 2,
    fixed_window = TRUE
  )
  train_bounds <- vapply(
    folds,
    function(fold) fold$train_idx,
    numeric(2)
  )
  expected <- cbind(
    c(1, 5),
    c(3, 7),
    c(5, 9),
    c(7, 11)
  )

  expect_identical(train_bounds, expected)
  for (fold in folds) {
    expect_length(fold$train_ts, 5)
  }
})


test_that("allow_partial retains the exact final incomplete horizon", {
  temp_data <- ts(seq_len(13), start = c(2000, 1), frequency = 4)

  complete <- ts_train_test_split(
    temp_data,
    initial = 6,
    horizon = 4,
    step = 4,
    allow_partial = FALSE
  )
  partial <- ts_train_test_split(
    temp_data,
    initial = 6,
    horizon = 4,
    step = 4,
    allow_partial = TRUE
  )
  final_fold <- partial[[2]]

  expect_length(complete, 1)
  expect_length(partial, 2)
  expect_identical(final_fold$train_idx, c(1, 10))
  expect_identical(final_fold$test_idx, c(11, 13))
  expect_identical(as.numeric(final_fold$test_ts), as.numeric(11:13))
  expect_identical(
    final_fold$test_range,
    stats::time(temp_data)[c(11, 13)]
  )
})


test_that("no complete horizon returns an empty fold list", {
  temp_data <- ts(seq_len(10), start = c(2000, 1), frequency = 4)

  folds <- ts_train_test_split(
    temp_data,
    initial = 8,
    horizon = 4,
    step = 1,
    allow_partial = FALSE
  )

  expect_identical(folds, list())
})


test_that("splits preserve frequency and multivariate exogenous values", {
  temp_data <- ts(seq_len(48), start = c(2000, 1), frequency = 24)
  exo_data <- ts(
    cbind(
      wind = 101:148,
      pressure = 201:248
    ),
    start = c(2000, 1),
    frequency = 24
  )

  folds <- ts_train_test_split(
    temp_data,
    exo_data = exo_data,
    initial = 24,
    horizon = 6,
    step = 12
  )
  first_fold <- folds[[1]]

  expect_length(folds, 2)
  expect_identical(stats::frequency(first_fold$train_ts), 24)
  expect_identical(stats::frequency(first_fold$test_ts), 24)
  expect_identical(stats::frequency(first_fold$exo_train_ts), 24)
  expect_identical(stats::frequency(first_fold$exo_test_ts), 24)
  expect_identical(
    as.numeric(first_fold$exo_train_ts[, "wind"]),
    as.numeric(101:124)
  )
  expect_identical(
    as.numeric(first_fold$exo_test_ts[, "pressure"]),
    as.numeric(225:230)
  )
  expect_identical(colnames(first_fold$exo_train_ts), c("wind", "pressure"))
  expect_identical(
    first_fold$train_range,
    stats::time(temp_data)[c(1, 24)]
  )
  expect_identical(
    first_fold$test_range,
    stats::time(temp_data)[c(25, 30)]
  )
})
