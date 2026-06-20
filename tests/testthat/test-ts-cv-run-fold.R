# test-ts-cv-run-fold.R

test_that("ts_cv_run_fold works without exogenous variables", {
  folds <- ts_train_test_split(
    temp_data = temp_ts_test,
    initial   = 60,
    horizon   = 12,
    step      = 12
  )

  res <- ts_cv_run_fold(folds[[1]])

  expect_type(res, "list")
  expect_true(res$converged)

  expect_s3_class(res$y_train, "ts")
  expect_s3_class(res$y_test, "ts")

  expect_false(is.null(res$y_pred))
})


test_that("ts_cv_run_fold errors on invalid fold", {
  expect_error(
    ts_cv_run_fold(NULL),
    "fold.*valid"
  )

  expect_error(
    ts_cv_run_fold(list()),
    "fold.*valid"
  )
})


test_that("ts_cv_run_fold validates model-control argument types", {
  folds <- ts_train_test_split(
    temp_data = temp_ts_test,
    initial   = 60,
    horizon   = 12,
    step      = 12
  )

  expect_error(
    ts_cv_run_fold(folds[[1]], ar_order = "1"),
    "ar_order.*numeric"
  )

  expect_error(
    ts_cv_run_fold(folds[[1]], use_season = 1),
    "use_season.*logical"
  )
})


test_that("ts_cv_run_fold rejects multivariate response series", {
  folds <- ts_train_test_split(
    temp_data = temp_ts_test,
    initial   = 60,
    horizon   = 12,
    step      = 12
  )

  bad_fold <- folds[[1]]
  bad_fold$train_ts <- ts(
    matrix(rnorm(length(bad_fold$train_ts) * 2), ncol = 2),
    start = start(bad_fold$train_ts),
    frequency = frequency(bad_fold$train_ts)
  )

  expect_error(
    ts_cv_run_fold(bad_fold),
    "fold\\$train_ts.*univariate"
  )
})
