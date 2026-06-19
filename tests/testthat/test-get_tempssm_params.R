# tests/testthat/test-get_tempssm_params.R

test_that("get_tempssm_params returns a list with expected names", {
  params <- get_tempssm_params(res_tempssm)

  expect_type(params, "list")

  expect_named(
    params,
    c("H", "Q_trend", "Q_season", "Q_ar", "ARs")
  )
})


test_that("get_tempssm_params returns numeric values", {
  params <- get_tempssm_params(res_tempssm)

  expect_type(params$H, "double")
  expect_type(params$Q_trend, "double")
  expect_type(params$Q_ar, "double")
  expect_type(params$ARs, "double")
})


test_that("Q_season is NA when seasonal is not used", {
  bad_res <- res_tempssm
  bad_res$use_season <- FALSE

  params <- get_tempssm_params(bad_res)

  expect_true(is.na(params$Q_season))
})


test_that("Q_season is computed when seasonal is used", {
  params <- get_tempssm_params(res_tempssm)

  expect_type(params$Q_season, "double")
  expect_false(is.na(params$Q_season))
})


test_that("parameters are exponentiated correctly", {
  pars <- res_tempssm$fit$optim.out$par

  ar_order <- res_tempssm$ar_order
  use_season <- res_tempssm$use_season

  if (use_season) {
    H_idx <- 4 + ar_order
  } else {
    H_idx <- 3 + ar_order
  }

  expected_H <- exp(pars[H_idx])

  params <- get_tempssm_params(res_tempssm)

  expect_identical(params$H, expected_H)
})


test_that("AR coefficients are transformed correctly", {
  params <- get_tempssm_params(res_tempssm)

  expect_length(params$ARs, res_tempssm$ar_order)
})


test_that("AR stationarity checker detects valid coefficients", {
  expect_true(.tempssm_is_stationary_ar(numeric(0)))
  expect_true(.tempssm_is_stationary_ar(0.5))
  expect_false(.tempssm_is_stationary_ar(1.1))
  expect_false(.tempssm_is_stationary_ar(NA_real_))
})


test_that("AR parameter transformation returns stationary coefficients", {
  ar_coefs <- .tempssm_transform_ar(c(10, -10, 2))

  expect_true(.tempssm_is_stationary_ar(ar_coefs))
})


test_that("get_tempssm_params checks input", {
  expect_error(
    get_tempssm_params(NULL),
    "`res` must be an object of class"
  )
})
