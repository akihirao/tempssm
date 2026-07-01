# tests/testthat/test-ts_CV_pipeline.R

test_that("lightweight CV pipeline works", {
  folds <- ts_train_test_split(
    temp_data = temp_ts_small,
    initial = 24,
    horizon = 12,
    step = 24
  )

  expect_length(folds, 1)

  cv_results <- ts_cv_run(
    folds,
    use_season = FALSE,
    parallel = FALSE,
    progress = FALSE
  )

  metrics <- lapply(cv_results, compute_cv_metrics)

  summary <- ts_cv_collect(cv_results, metrics)

  expect_s3_class(summary, "tbl_df")
  expect_identical(nrow(summary), 1L)
})


test_that("ts_cv_run validates input types before execution", {
  expect_error(
    ts_cv_run(1),
    "folds.*list"
  )

  expect_error(
    ts_cv_run(list(), ar_order = "1"),
    "ar_order.*numeric"
  )

  expect_error(
    ts_cv_run(list(), use_season = 1),
    "use_season.*logical"
  )

  expect_error(
    ts_cv_run(list(), parallel = 1),
    "parallel.*logical"
  )

  expect_error(
    ts_cv_run(list(), workers = "2"),
    "workers.*numeric"
  )

  expect_error(
    ts_cv_run(list(), progress = 1),
    "progress.*logical"
  )
})


test_that("ts_cv_run preserves fold order and forwards model controls", {
  folds <- rev(ts_train_test_split(
    temp_data = temp_ts_small,
    initial = 24,
    horizon = 12,
    step = 12
  ))

  testthat::local_mocked_bindings(
    ts_cv_run_fold = function(fold, ar_order, use_season, marginal) {
      list(
        fold = fold$fold,
        ar_order = ar_order,
        use_season = use_season,
        marginal = marginal
      )
    },
    .package = "tempssm"
  )

  res <- ts_cv_run(
    folds,
    ar_order = 2,
    use_season = FALSE,
    marginal = TRUE,
    parallel = FALSE,
    workers = 1,
    progress = FALSE
  )

  expect_identical(vapply(res, `[[`, numeric(1), "fold"), c(2, 1))
  expect_identical(vapply(res, `[[`, numeric(1), "ar_order"), c(2, 2))
  expect_false(any(vapply(res, `[[`, logical(1), "use_season")))
  expect_true(all(vapply(res, `[[`, logical(1), "marginal")))
})


test_that("ts_cv_run accepts an empty fold list", {
  res <- ts_cv_run(
    list(),
    parallel = FALSE,
    workers = 1,
    progress = FALSE
  )

  expect_identical(res, list())
})


test_that("ts_cv_run restores the previous future plan", {
  original_plan <- future::plan()
  withr::defer(future::plan(original_plan))
  future::plan(future::multisession, workers = 1)
  expected_plan <- future::plan()

  ts_cv_run(
    list(),
    parallel = FALSE,
    workers = 1,
    progress = FALSE
  )

  expect_identical(future::plan(), expected_plan)
})


test_that("temporary CV plans are restored after errors", {
  original_plan <- future::plan()
  withr::defer(future::plan(original_plan))
  future::plan(future::multisession, workers = 1)
  expected_plan <- future::plan()

  expect_error(
    .with_ts_cv_plan(
      parallel = FALSE,
      workers = 1,
      code = function() stop("synthetic execution failure")
    ),
    "synthetic execution failure"
  )

  expect_identical(future::plan(), expected_plan)
})


test_that("ts_cv_run validates folds before execution", {
  bad_fold <- list(
    fold = 1L,
    train_ts = ts(1:12, frequency = 12),
    test_ts = ts(13:18, start = c(2, 1), frequency = 12),
    exo_train_ts = ts(1:12, frequency = 12),
    exo_test_ts = NULL
  )

  expect_error(
    ts_cv_run(
      list(bad_fold),
      parallel = FALSE,
      workers = 1
    ),
    "Exogenous fold series.*both"
  )
})


test_that("ts_cv_run supports progress reporting", {
  testthat::skip_if_not_installed("progressr")
  folds <- ts_train_test_split(
    temp_data = temp_ts_small,
    initial = 24,
    horizon = 12,
    step = 24
  )

  testthat::local_mocked_bindings(
    ts_cv_run_fold = function(fold, ar_order, use_season, marginal) {
      list(fold = fold$fold, converged = TRUE)
    },
    .package = "tempssm"
  )

  res <- ts_cv_run(
    folds,
    parallel = FALSE,
    workers = 1,
    progress = TRUE
  )

  expect_length(res, 1L)
  expect_true(res[[1]]$converged)
})
