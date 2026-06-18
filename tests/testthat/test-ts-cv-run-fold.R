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
    ts_cv_run_fold(NULL)
  )

  expect_error(
    ts_cv_run_fold(list())
  )
})
