# tests/testthat/test-ts_CV_pipeline.R

test_that("lightweight CV pipeline works", {
  folds <- ts_train_test_split(
    temp_data = temp_ts_small,
    initial = 24,
    horizon = 12,
    step = 24
  )

  expect_equal(length(folds), 1)

  cv_results <- ts_cv_run(
    folds,
    use_season = FALSE,
    parallel = FALSE,
    progress = FALSE
  )

  metrics <- lapply(cv_results, compute_cv_metrics)

  summary <- ts_cv_collect(cv_results, metrics)

  expect_s3_class(summary, "tbl_df")
  expect_equal(nrow(summary), 1)
})
