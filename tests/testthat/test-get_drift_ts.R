# tests/testthat/test-get_drift_ts.R

test_that("get_drift_ts basic structure", {
  ts_obj <- get_drift_ts(res_tempssm)

  check_ts_basic(ts_obj, res_tempssm)
  check_ts_univariate(ts_obj)
})


test_that("get_drift_ts CI structure", {
  ts_ci <- get_drift_ts(res_tempssm, ci = TRUE)

  check_ts_ci_structure(ts_ci, "drift")
})

test_that("get_drift_ts CI values scaled", {
  ts_ci <- get_drift_ts(res_tempssm, ci = TRUE)

  ci_obj <- stats::confint(res_tempssm$kfs)
  freq <- frequency(res_tempssm$temp_data)

  check_ts_ci_values(ts_ci, ci_obj, "slope", scale = freq)
})


test_that("get_drift_ts preserves time attributes", {
  ts_obj <- get_drift_ts(res_tempssm)

  expect_equal(start(ts_obj), start(res_tempssm$temp_data))
  expect_equal(frequency(ts_obj), frequency(res_tempssm$temp_data))
})


test_that("CI is scaled by frequency", {

  freq <- frequency(res_tempssm$temp_data)

  ts_ci <- get_drift_ts(res_tempssm, ci = TRUE)

  ci_obj <- stats::confint(res_tempssm$kfs)

  expect_equal(
    ts_ci[, "lwr"],
    ci_obj$slope[, "lwr"] * freq
  )

  expect_equal(
    ts_ci[, "upr"],
    ci_obj$slope[, "upr"] * freq
  )
})


test_that("CI output has correct structure", {
  ts_ci <- get_drift_ts(res_tempssm, ci = TRUE)

  expect_equal(NCOL(ts_ci), 3)
  expect_named(as.data.frame(ts_ci), c("drift", "lwr", "upr"))
})


test_that("errors when slope is missing in alphahat", {

  bad_res <- res_tempssm
  bad_res$kfs$alphahat <- matrix(rnorm(10), ncol = 1)
  colnames(bad_res$kfs$alphahat) <- "level"

  expect_error(
    get_drift_ts(bad_res),
    "Drift \\(slope\\) component not found"
  )
})



test_that("errors when smoothing results are missing", {

  bad_res <- res_tempssm
  bad_res$kfs$alphahat <- NULL

  expect_error(
    get_drift_ts(bad_res),
    "Drift \\(slope\\) component not found"
  )
})



test_that("errors when slope missing in confidence intervals", {

  bad_res <- res_tempssm

  testthat::local_mocked_bindings(
    confint = function(...) list(level = matrix(0)),
    .package = "stats"
  )

  expect_error(
    get_drift_ts(bad_res, ci = TRUE),
    "Slope component not found in confidence intervals"
  )
})



test_that("ci = FALSE returns univariate ts", {
  ts_obj <- get_drift_ts(res_tempssm, ci = FALSE)

  expect_equal(NCOL(ts_obj), 1)
})


test_that("ci_level boundary values", {
  expect_error(get_drift_ts(res_tempssm, ci = TRUE, ci_level = 0))
  expect_error(get_drift_ts(res_tempssm, ci = TRUE, ci_level = 1))
})


