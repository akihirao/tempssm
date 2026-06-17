# tests/testthat/test-get_params.R

test_that("get_params returns a list with expected names", {
  params <- get_params(res_tempssm)

  expect_type(params, "list")

  expect_named(
    params,
    c("H", "Q_trend", "Q_season", "Q_ar", "ARs")
  )
})


test_that("get_params returns numeric values", {
  params <- get_params(res_tempssm)

  expect_true(is.numeric(params$H))
  expect_true(is.numeric(params$Q_trend))
  expect_true(is.numeric(params$Q_ar))
  expect_true(is.numeric(params$ARs))
})


test_that("Q_season is NA when seasonal is not used", {
  bad_res <- res_tempssm
  bad_res$use_season <- FALSE

  params <- get_params(bad_res)

  expect_true(is.na(params$Q_season))
})


test_that("Q_season is computed when seasonal is used", {
  params <- get_params(res_tempssm)

  expect_true(is.numeric(params$Q_season))
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

  params <- get_params(res_tempssm)

  expect_equal(params$H, expected_H)
})


test_that("AR coefficients are transformed correctly", {
  params <- get_params(res_tempssm)

  expect_equal(length(params$ARs), res_tempssm$ar_order)
})


test_that("get_params checks input", {
  expect_error(
    get_params(NULL),
    "`res` must be an object of class"
  )
})
