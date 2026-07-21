# tests/testthat/test-validate_tempssm_for_logLik.R

test_that(".validate_tempssm_for_logLik works for valid object", {
  expect_invisible(
    out <- .validate_tempssm_for_logLik(res_tempssm)
  )

  expect_identical(out, res_tempssm)
})


test_that(".validate_tempssm_for_logLik errors for invalid class", {
  expect_error(
    .validate_tempssm_for_logLik(NULL),
    "`res` must be an object of class"
  )
})


test_that(".validate_tempssm_for_logLik errors when not converged", {
  bad_res <- res_tempssm
  bad_res$converged <- FALSE

  expect_error(
    .validate_tempssm_for_logLik(bad_res),
    "did not converge"
  )
})


test_that(".validate_tempssm_for_logLik errors when model is missing", {
  bad_res <- res_tempssm
  bad_res$model <- NULL

  expect_error(
    .validate_tempssm_for_logLik(bad_res),
    "fitted model is missing"
  )
})


test_that(".validate_tempssm_for_logLik errors when fit is missing", {
  bad_res <- res_tempssm
  bad_res$fit <- NULL

  expect_error(
    .validate_tempssm_for_logLik(bad_res),
    "fitted model is missing"
  )
})


test_that(".validate_tempssm_for_logLik errors if model or fit is NULL", {
  bad_res1 <- res_tempssm
  bad_res1$model <- NULL

  bad_res2 <- res_tempssm
  bad_res2$fit <- NULL

  expect_error(
    .validate_tempssm_for_logLik(bad_res1),
    "fitted model is missing"
  )
  expect_error(
    .validate_tempssm_for_logLik(bad_res2),
    "fitted model is missing"
  )
})


test_that(".validate_tempssm_for_logLik does not modify object", {
  res_copy <- res_tempssm

  expect_error(
    tryCatch(
      .validate_tempssm_for_logLik(res_copy),
      error = function(e) NULL
    ),
    NA
  )

  expect_identical(res_copy, res_tempssm)
})


test_that(".validate_tempssm_for_logLik requires TRUE convergence", {
  bad_res <- res_tempssm
  bad_res$converged <- NA # subtle case

  expect_error(
    .validate_tempssm_for_logLik(bad_res),
    "did not converge"
  )
})
