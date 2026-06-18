# tests/testthat/test-robustness-non-convergence.R
# Test robustness of functions when tempssm models do not converge

test_that("get_exo_coef returns NULL for non-converged models", {
  res <- res_tempssm_exo
  res$converged <- FALSE

  result <- get_exo_coef(res)
  expect_null(result)
})


test_that("get_tempssm_params handles missing fit results gracefully", {
  res <- res_tempssm
  res$fit <- NULL

  expect_error(
    get_tempssm_params(res),
    "Fitted parameters are not available"
  )
})


test_that("get_tempssm_params handles missing optim.out gracefully", {
  res <- res_tempssm
  res$fit$optim.out <- NULL

  expect_error(
    get_tempssm_params(res),
    "Fitted parameters are not available"
  )
})


test_that("summary.tempssm handles non-converged models", {
  res <- res_tempssm
  res$converged <- FALSE
  res$fit <- NULL

  expect_error(
    summary(res),
    "fitted model results are missing"
  )
})


test_that("summary.tempssm handles missing fit results", {
  res <- res_tempssm
  res$fit <- NULL

  expect_error(
    summary(res),
    "fitted model results are missing"
  )
})


test_that("summary.tempssm handles missing optim.out gracefully", {
  res <- res_tempssm
  res$fit$optim.out <- NULL

  expect_error(
    summary(res),
    "fitted model results are missing"
  )
})


test_that("get_level_ts works even if model did not converge", {
  # Note: get_level_ts only checks alphahat exists,
  # it does not care about convergence status
  res <- res_tempssm
  res$converged <- FALSE

  ts_obj <- get_level_ts(res)

  expect_s3_class(ts_obj, "ts")
  expect_identical(NCOL(ts_obj), 1L)
})


test_that("get_drift_ts works even if model did not converge", {
  res <- res_tempssm
  res$converged <- FALSE

  ts_obj <- get_drift_ts(res)

  expect_s3_class(ts_obj, "ts")
  expect_identical(NCOL(ts_obj), 1L)
})


test_that("get_season_ts returns NULL for non-seasonal converged model", {
  res <- res_tempssm
  res$use_season <- FALSE

  expect_error(
    get_season_ts(res),
    "Seasonal component is not included"
  )
})


test_that("logLik.tempssm handles non-converged models", {
  res <- res_tempssm
  res$converged <- FALSE

  expect_error(
    logLik(res),
    "model did not converge"
  )
})


test_that("AIC.tempssm handles non-converged models", {
  res <- res_tempssm
  res$converged <- FALSE

  expect_error(
    AIC(res),
    "model did not converge"
  )
})


test_that("logLik.tempssm handles completely missing fitted model", {
  res <- res_tempssm
  res$model <- NULL
  res$fit <- NULL

  expect_error(
    logLik(res),
    "fitted model is missing"
  )
})
