# tests/testthat/test-ts_slice.R

test_that("ts_train_test_split returns list of folds", {
  folds <- ts_train_test_split(temp_ts_test)

  expect_type(folds, "list")
  expect_true(length(folds) > 0)

  expect_true(all(sapply(folds, function(x) is.list(x))))
})


test_that("number of folds matches expected", {
  n <- length(temp_ts_test)
  initial <- 60
  horizon <- 12
  step <- 12

  expected <- floor((n - initial) / step)

  folds <- ts_train_test_split(
    temp_ts_test,
    initial = initial,
    horizon = horizon,
    step = step
  )

  expect_equal(length(folds), expected)
})


test_that("each fold has required elements", {
  folds <- ts_train_test_split(temp_ts_test)

  f1 <- folds[[1]]

  expect_named(
    f1,
    c(
      "fold", "train_ts", "test_ts",
      "exo_train_ts", "exo_test_ts",
      "train_idx", "test_idx",
      "train_range", "test_range"
    )
  )
})


test_that("train and test indices are consistent", {
  folds <- ts_train_test_split(temp_ts_test)

  f1 <- folds[[1]]

  t_idx <- f1$train_idx
  s_idx <- f1$test_idx

  expect_true(t_idx[1] < t_idx[2])
  expect_true(s_idx[1] < s_idx[2])

  expect_equal(s_idx[1], t_idx[2] + 1)
})


test_that("train and test lengths match indices", {
  folds <- ts_train_test_split(temp_ts_test)

  f <- folds[[1]]

  expect_equal(
    length(f$train_ts),
    diff(f$train_idx) + 1
  )

  expect_equal(
    length(f$test_ts),
    diff(f$test_idx) + 1
  )
})


test_that("fixed vs expanding window works correctly", {
  f_exp <- ts_train_test_split(temp_ts_test,
    fixed_window = FALSE
  )

  f_fix <- ts_train_test_split(temp_ts_test,
    fixed_window = TRUE
  )

  # expanding → 長くなる
  expect_true(
    length(f_exp[[2]]$train_ts) >
      length(f_exp[[1]]$train_ts)
  )

  # fixed → 一定
  expect_equal(
    length(f_fix[[2]]$train_ts),
    length(f_fix[[1]]$train_ts)
  )
})


test_that("allow_partial changes final fold", {
  f1 <- ts_train_test_split(
    temp_ts_test,
    allow_partial = FALSE
  )

  f2 <- ts_train_test_split(
    temp_ts_test,
    allow_partial = TRUE
  )

  expect_true(length(f2) >= length(f1))
})


test_that("works with exogenous variables", {
  folds <- ts_train_test_split(
    temp_ts_test,
    exo_data = exo_ts_test
  )

  f1 <- folds[[1]]

  expect_s3_class(f1$exo_train_ts, "ts")
  expect_s3_class(f1$exo_test_ts, "ts")

  expect_equal(
    length(f1$exo_train_ts),
    length(f1$train_ts)
  )
})


test_that("exo_data is NULL when not provided", {
  folds <- ts_train_test_split(temp_ts_test)

  expect_null(folds[[1]]$exo_train_ts)
  expect_null(folds[[1]]$exo_test_ts)
})
