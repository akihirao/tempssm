test_that(".build_ts_split_fold constructs every fold field", {
  temp_data <- ts(seq_len(12), start = c(2000, 1), frequency = 4)
  exo_data <- ts(
    cbind(
      wind = 101:112,
      pressure = 201:212
    ),
    start = c(2000, 1),
    frequency = 4
  )

  fold <- .build_ts_split_fold(
    fold_id = 2,
    temp_data = temp_data,
    exo_data = exo_data,
    train_start = 3,
    train_end = 8,
    test_start = 9,
    test_end = 11
  )

  expect_named(
    fold,
    c(
      "fold",
      "train_ts",
      "test_ts",
      "exo_train_ts",
      "exo_test_ts",
      "train_idx",
      "test_idx",
      "train_range",
      "test_range"
    )
  )
  expect_identical(fold$fold, 2)
  expect_identical(fold$train_idx, c(3, 8))
  expect_identical(fold$test_idx, c(9, 11))
  expect_identical(as.numeric(fold$train_ts), as.numeric(3:8))
  expect_identical(as.numeric(fold$test_ts), as.numeric(9:11))
  expect_identical(
    as.numeric(fold$exo_train_ts[, "wind"]),
    as.numeric(103:108)
  )
  expect_identical(
    as.numeric(fold$exo_test_ts[, "pressure"]),
    as.numeric(209:211)
  )
  expect_identical(fold$train_range, stats::time(temp_data)[c(3, 8)])
  expect_identical(fold$test_range, stats::time(temp_data)[c(9, 11)])
})


test_that(".build_ts_split_fold preserves absent exogenous data", {
  temp_data <- ts(seq_len(8), frequency = 4)

  fold <- .build_ts_split_fold(
    fold_id = 1,
    temp_data = temp_data,
    exo_data = NULL,
    train_start = 1,
    train_end = 4,
    test_start = 5,
    test_end = 8
  )

  expect_null(fold$exo_train_ts)
  expect_null(fold$exo_test_ts)
})
