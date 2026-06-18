# tests/testthat/test-get_level_ts.R

test_that("get_level_ts returns a ts object", {
  ts_obj <- get_level_ts(res_tempssm)

  expect_s3_class(ts_obj, "ts")
  expect_identical(NCOL(ts_obj), 1L)
})


test_that("get_level_ts returns CI columns when ci = TRUE", {
  ts_ci <- get_level_ts(res_tempssm, ci = TRUE)

  expect_s3_class(ts_ci, "ts")
  expect_identical(colnames(ts_ci), c("level", "lwr", "upr"))
})


test_that("get_level_ts checks inputs correctly", {
  expect_error(
    get_level_ts(NULL),
    "`res` must be an object of class 'tempssm'."
  )

  expect_error(
    get_level_ts(res_tempssm, ci = TRUE, ci_level = 1.2),
    "`ci_level` must be a numeric value between 0 and 1."
  )
})


test_that("get_level_ts preserves time attributes", {
  ts_obj <- get_level_ts(res_tempssm)

  expect_identical(start(ts_obj), start(res_tempssm$temp_data))
  expect_identical(frequency(ts_obj), frequency(res_tempssm$temp_data))
})


test_that("get_level_ts returns correct number of columns with CI", {
  ts_ci <- get_level_ts(res_tempssm, ci = TRUE)

  expect_identical(NCOL(ts_ci), 3L)
})


test_that("CI output preserves time attributes", {
  ts_ci <- get_level_ts(res_tempssm, ci = TRUE)

  expect_identical(start(ts_ci), start(res_tempssm$temp_data))
  expect_identical(frequency(ts_ci), frequency(res_tempssm$temp_data))
})


test_that("errors when level component is missing", {
  bad_res <- res_tempssm
  bad_res$kfs$alphahat <- NULL

  expect_error(
    get_level_ts(bad_res),
    "Level component not found"
  )
})


test_that("errors when level column is missing in alphahat", {
  bad_res <- res_tempssm
  bad_res$kfs$alphahat <- matrix(rnorm(10), ncol = 1)
  colnames(bad_res$kfs$alphahat) <- "trend" # no existence of "level"
  expect_error(
    get_level_ts(bad_res),
    "Level component not found"
  )
})


test_that("errors when level missing in confidence intervals", {
  bad_res <- res_tempssm

  original_confint <- stats::confint

  testthat::local_mocked_bindings(
    confint = function(...) list(other = matrix(0)),
    .package = "stats"
  )

  expect_error(
    get_level_ts(bad_res, ci = TRUE),
    "Level component not found in confidence intervals"
  )
})


test_that("ci = FALSE returns univariate ts", {
  ts_obj <- get_level_ts(res_tempssm, ci = FALSE)

  expect_identical(NCOL(ts_obj), 1L)
})


test_that("ci = FALSE returns univariate ts", {
  ts_obj <- get_level_ts(res_tempssm, ci = FALSE)

  expect_identical(NCOL(ts_obj), 1L)
})


test_that("ci_level boundary values", {
  expect_error(get_level_ts(res_tempssm, ci = TRUE, ci_level = 0))
  expect_error(get_level_ts(res_tempssm, ci = TRUE, ci_level = 1))
})
