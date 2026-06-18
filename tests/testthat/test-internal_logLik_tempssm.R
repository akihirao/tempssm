# tests/testthat/test-internal_logLik_tempssm.R

test_that(".internal_logLik_tempssm returns correct structure", {
  res <- tempssm:::.internal_logLik_tempssm(res_tempssm)

  expect_type(res, "list")

  expect_named(
    res,
    c("logLik", "df", "nobs")
  )
})


test_that(".internal_logLik_tempssm returns numeric values", {
  res <- tempssm:::.internal_logLik_tempssm(res_tempssm)

  expect_type(res$logLik, "double")
  expect_type(res$df, "integer")
  expect_type(res$nobs, "integer")
})


test_that("nobs equals length of temp_data", {
  out <- tempssm:::.internal_logLik_tempssm(res_tempssm)

  expect_identical(out$nobs, length(res_tempssm$temp_data))
})


test_that("df equals number of parameters without exogenous variables", {
  out <- tempssm:::.internal_logLik_tempssm(res_tempssm)

  expected_df <- length(res_tempssm$fit$optim.out$par)

  expect_identical(out$df, expected_df)
})


test_that("df includes exogenous variables when present", {
  res <- res_tempssm_exo

  out <- tempssm:::.internal_logLik_tempssm(res)

  expected_df <- length(res$fit$optim.out$par) + ncol(res$exogenous_data)

  expect_identical(out$df, expected_df)
})


test_that("logLik matches stats::logLik output", {
  out <- tempssm:::.internal_logLik_tempssm(res_tempssm)

  expected <- as.numeric(stats::logLik(res_tempssm$model))

  expect_identical(out$logLik, expected)
})


test_that(".internal_logLik_tempssm errors for invalid input", {
  expect_error(
    tempssm:::.internal_logLik_tempssm(NULL),
    "`res` must be an object of class"
  )
})


test_that(".internal_logLik_tempssm errors when not converged", {
  bad_res <- res_tempssm
  bad_res$converged <- FALSE

  expect_error(
    tempssm:::.internal_logLik_tempssm(bad_res),
    "did not converge"
  )
})


test_that("errors when logLik extraction fails", {
  bad_res <- res_tempssm

  # break model
  bad_res$model <- list()

  expect_error(
    tempssm:::.internal_logLik_tempssm(bad_res),
    "Failed to extract log-likelihood"
  )
})


test_that("handles exogenous with zero columns", {
  res <- res_tempssm
  res$exogenous_data <- matrix(nrow = length(res$temp_data), ncol = 0)

  out <- tempssm:::.internal_logLik_tempssm(res)

  expect_identical(
    out$df,
    length(res$fit$optim.out$par)
  )
})


test_that("function does not modify input object", {
  res_copy <- res_tempssm

  out <- tempssm:::.internal_logLik_tempssm(res_copy)

  expect_identical(res_copy, res_tempssm)
})
