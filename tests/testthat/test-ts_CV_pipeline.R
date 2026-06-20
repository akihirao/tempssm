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
