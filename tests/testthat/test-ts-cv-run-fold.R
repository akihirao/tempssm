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


test_that("ts_cv_run_fold returns the documented success structure", {
  fold <- list(
    fold = 3L,
    train_ts = ts(1:24, start = c(2000, 1), frequency = 12),
    test_ts = ts(25:30, start = c(2002, 1), frequency = 12),
    exo_train_ts = NULL,
    exo_test_ts = NULL
  )
  fitted_model <- structure(list(id = "fitted"), class = "mock_model")
  fitted_result <- structure(
    list(converged = TRUE, model = fitted_model),
    class = "tempssm"
  )

  testthat::local_mocked_bindings(
    .fit_tempssm_safe = function(...) fitted_result,
    .predict_no_exo = function(model, h) seq_len(h),
    .package = "tempssm"
  )

  res <- ts_cv_run_fold(fold, ar_order = 2, use_season = FALSE)

  expect_named(
    res,
    c("fold", "converged", "y_train", "y_test", "y_pred", "model")
  )
  expect_identical(res$fold, 3L)
  expect_true(res$converged)
  expect_identical(res$y_train, fold$train_ts)
  expect_identical(res$y_test, fold$test_ts)
  expect_identical(res$y_pred, seq_len(NROW(fold$test_ts)))
  expect_identical(res$model, fitted_model)
})


test_that("ts_cv_run_fold dispatches exogenous forecasts", {
  fold <- list(
    fold = 4L,
    train_ts = ts(1:24, start = c(2000, 1), frequency = 12),
    test_ts = ts(25:30, start = c(2002, 1), frequency = 12),
    exo_train_ts = ts(
      matrix(1:24, ncol = 1, dimnames = list(NULL, "x")),
      start = c(2000, 1),
      frequency = 12
    ),
    exo_test_ts = ts(
      matrix(25:30, ncol = 1, dimnames = list(NULL, "x")),
      start = c(2002, 1),
      frequency = 12
    )
  )
  fitted_result <- structure(
    list(converged = TRUE, model = list(id = "fitted")),
    class = "tempssm"
  )
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    .fit_tempssm_safe = function(...) fitted_result,
    .predict_with_exo = function(res, y_train_named, y_test_named,
                                 exo_test, ar_order, use_season) {
      observed$exo_test <- exo_test
      observed$ar_order <- ar_order
      observed$use_season <- use_season
      rep(10, NROW(y_test_named))
    },
    .package = "tempssm"
  )

  res <- ts_cv_run_fold(fold, ar_order = 2, use_season = FALSE)

  expect_identical(res$y_pred, rep(10, NROW(fold$test_ts)))
  expect_identical(observed$exo_test, fold$exo_test_ts)
  expect_identical(observed$ar_order, 2)
  expect_false(observed$use_season)
})


test_that("ts_cv_run_fold forwards marginal likelihood control", {
  fold <- list(
    fold = 6L,
    train_ts = ts(1:24, frequency = 12),
    test_ts = ts(25:30, start = c(3, 1), frequency = 12),
    exo_train_ts = NULL,
    exo_test_ts = NULL
  )
  observed <- new.env(parent = emptyenv())
  fitted_result <- structure(
    list(converged = TRUE, model = list(id = "fitted")),
    class = "tempssm"
  )

  testthat::local_mocked_bindings(
    .fit_tempssm_safe = function(...) {
      observed$args <- list(...)
      fitted_result
    },
    .predict_no_exo = function(model, h) seq_len(h),
    .package = "tempssm"
  )

  ts_cv_run_fold(fold, marginal = TRUE)

  expect_true(observed$args$marginal)
})


test_that("ts_cv_run_fold returns the documented failure structure", {
  fold <- list(
    fold = 5L,
    train_ts = ts(1:24, frequency = 12),
    test_ts = ts(25:30, start = c(3, 1), frequency = 12),
    exo_train_ts = NULL,
    exo_test_ts = NULL
  )

  testthat::local_mocked_bindings(
    .fit_tempssm_safe = function(...) NULL,
    .package = "tempssm"
  )

  expect_warning(
    res <- ts_cv_run_fold(fold),
    "Model did not converge for fold 5"
  )

  expect_named(
    res,
    c("fold", "converged", "y_train", "y_test", "y_pred", "model")
  )
  expect_identical(res$fold, 5L)
  expect_false(res$converged)
  expect_identical(res$y_train, fold$train_ts)
  expect_identical(res$y_test, fold$test_ts)
  expect_null(res$y_pred)
  expect_null(res$model)
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

  invalid_id <- list(
    fold = 0,
    train_ts = ts(1:12, frequency = 12),
    test_ts = ts(13:18, start = c(2, 1), frequency = 12)
  )
  expect_error(
    ts_cv_run_fold(invalid_id),
    "fold\\$fold.*positive integer"
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


test_that("ts_cv_run_fold requires paired exogenous series", {
  folds <- ts_train_test_split(
    temp_data = temp_ts_test,
    exo_data = exo_ts_test,
    initial = 60,
    horizon = 12,
    step = 12
  )
  bad_fold <- folds[[1]]
  bad_fold$exo_test_ts <- NULL

  expect_error(
    ts_cv_run_fold(bad_fold),
    "Exogenous fold series.*both"
  )
})


test_that("ts_cv_run_fold validates test exogenous alignment", {
  folds <- ts_train_test_split(
    temp_data = temp_ts_test,
    exo_data = exo_ts_test,
    initial = 60,
    horizon = 12,
    step = 12
  )
  bad_fold <- folds[[1]]
  bad_fold$exo_test_ts <- stats::window(
    bad_fold$exo_test_ts,
    end = stats::time(bad_fold$exo_test_ts)[6]
  )

  expect_error(
    ts_cv_run_fold(bad_fold),
    "Length of.*exo_data.*temp_data"
  )
})
