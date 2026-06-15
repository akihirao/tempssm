# tests/testthat/test-ts-cv-collect.R

test_that("ts_cv_collect returns a tibble", {

  cv_results <- list(
    list(fold = 1, converged = TRUE),
    list(fold = 2, converged = FALSE)
  )

  metrics <- list(
    list(MAE = 1, MASE_naive = 0.5, MASE_seasonal = 0.7),
    list(MAE = 2, MASE_naive = 0.6, MASE_seasonal = 0.8)
  )

  res <- ts_cv_collect(cv_results, metrics)

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 2)
})


test_that("ts_cv_collect returns correct column names", {

  cv_results <- list(
    list(fold = 1, converged = TRUE)
  )

  metrics <- list(
    list(MAE = 1, MASE_naive = 0.5, MASE_seasonal = 0.7)
  )

  res <- ts_cv_collect(cv_results, metrics)

  expect_named(
    res,
    c("fold", "converged", "MAE", "MASE_naive", "MASE_seasonal")
  )
})


test_that("values are correctly extracted", {

  cv_results <- list(
    list(fold = 1, converged = TRUE),
    list(fold = 2, converged = FALSE)
  )

  metrics <- list(
    list(MAE = 10, MASE_naive = 1.1, MASE_seasonal = 1.2),
    list(MAE = 20, MASE_naive = 2.1, MASE_seasonal = 2.2)
  )

  res <- ts_cv_collect(cv_results, metrics)

  expect_equal(res$fold, c(1, 2))
  expect_equal(res$converged, c(1, 0))  # as.numeric(TRUE/FALSE)
  expect_equal(res$MAE, c(10, 20))
})


test_that("outputs are numeric", {

  cv_results <- list(
    list(fold = 1, converged = TRUE)
  )

  metrics <- list(
    list(MAE = 1, MASE_naive = 0.5, MASE_seasonal = 0.7)
  )

  res <- ts_cv_collect(cv_results, metrics)

  expect_type(res$fold, "double")
  expect_type(res$converged, "double")
})


test_that("errors when inputs are not lists", {

  expect_error(
    ts_cv_collect(1, list()),
    "must both be lists"
  )

  expect_error(
    ts_cv_collect(list(), 1),
    "must both be lists"
  )
})


test_that("errors when lengths differ", {

  cv_results <- list(list(fold=1, converged=TRUE))
  metrics <- list()

  expect_error(
    ts_cv_collect(cv_results, metrics),
    "same length"
  )
})


test_that("handles NA values", {

  cv_results <- list(
    list(fold = 1, converged = TRUE)
  )

  metrics <- list(
    list(MAE = NA, MASE_naive = NA, MASE_seasonal = NA)
  )

  res <- ts_cv_collect(cv_results, metrics)

  expect_true(is.na(res$MAE))
})

