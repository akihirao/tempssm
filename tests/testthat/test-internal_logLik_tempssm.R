# tests/testthat/test-internal_logLik_tempssm.R

test_that(".internal_logLik_tempssm returns correct structure", {
  res <- .internal_logLik_tempssm(res_tempssm)

  expect_type(res, "list")

  expect_named(
    res,
    c("logLik", "df", "nobs", "marginal")
  )
})


test_that(".internal_logLik_tempssm returns numeric values", {
  res <- .internal_logLik_tempssm(res_tempssm)

  expect_type(res$logLik, "double")
  expect_type(res$df, "integer")
  expect_type(res$nobs, "integer")
  expect_type(res$marginal, "logical")
})


test_that("nobs equals length of temp_data", {
  out <- .internal_logLik_tempssm(res_tempssm)

  expect_identical(out$nobs, length(res_tempssm$temp_data))
})


test_that("df equals number of parameters without exogenous variables", {
  out <- .internal_logLik_tempssm(res_tempssm)

  expected_df <- length(res_tempssm$fit$optim.out$par)

  expect_identical(out$df, expected_df)
})


test_that("df includes exogenous variables when present", {
  res <- res_tempssm_exo

  out <- .internal_logLik_tempssm(res)

  expected_df <- length(res$fit$optim.out$par) + ncol(res$exogenous_data)

  expect_identical(out$df, expected_df)
})


test_that("logLik matches stats::logLik output", {
  out <- .internal_logLik_tempssm(res_tempssm)

  expected <- as.numeric(stats::logLik(
    res_tempssm$model,
    marginal = res_tempssm$marginal
  ))

  expect_identical(out$logLik, expected)
})


test_that("logLik can evaluate diffuse and marginal likelihoods", {
  diffuse <- .internal_logLik_tempssm(res_tempssm, marginal = FALSE)
  marginal <- .internal_logLik_tempssm(res_tempssm, marginal = TRUE)

  expect_identical(
    diffuse$logLik,
    as.numeric(stats::logLik(res_tempssm$model, marginal = FALSE))
  )
  expect_identical(
    marginal$logLik,
    as.numeric(stats::logLik(res_tempssm$model, marginal = TRUE))
  )
  expect_false(diffuse$marginal)
  expect_true(marginal$marginal)
})


test_that("logLik uses stored setting and supports legacy objects", {
  marginal_res <- res_tempssm
  marginal_res$marginal <- TRUE
  expect_true(attr(logLik(marginal_res), "marginal"))

  legacy_res <- res_tempssm
  legacy_res$marginal <- NULL
  expect_false(attr(logLik(legacy_res), "marginal"))
})


test_that("logLik validates explicit marginal setting", {
  expect_error(
    logLik(res_tempssm, marginal = NA),
    "marginal.*logical scalar"
  )
  expect_error(
    logLik(res_tempssm, marginal = c(TRUE, FALSE)),
    "marginal.*length one"
  )
})


test_that(".internal_logLik_tempssm errors for invalid input", {
  expect_error(
    .internal_logLik_tempssm(NULL),
    "`res` must be an object of class"
  )
})


test_that(".internal_logLik_tempssm errors when not converged", {
  bad_res <- res_tempssm
  bad_res$converged <- FALSE

  expect_error(
    .internal_logLik_tempssm(bad_res),
    "did not converge"
  )
})


test_that("errors when logLik extraction fails", {
  bad_res <- res_tempssm

  # break model
  bad_res$model <- list()

  expect_error(
    .internal_logLik_tempssm(bad_res),
    "Failed to extract log-likelihood"
  )
})


test_that("handles exogenous with zero columns", {
  res <- res_tempssm
  res$exogenous_data <- matrix(nrow = length(res$temp_data), ncol = 0)

  out <- .internal_logLik_tempssm(res)

  expect_identical(
    out$df,
    length(res$fit$optim.out$par)
  )
})


test_that("function does not modify input object", {
  res_copy <- res_tempssm

  out <- .internal_logLik_tempssm(res_copy)

  expect_identical(res_copy, res_tempssm)
})
