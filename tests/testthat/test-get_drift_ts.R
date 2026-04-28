#tests/testthat/test-get_drift_ts.R

test_that("get_drift_ts returns a ts object", {

  ts_obj <- get_drift_ts(res_tempssm)

  expect_s3_class(ts_obj, "ts")
  expect_equal(NCOL(ts_obj), 1)
})


test_that("get_drift_ts returns CI columns when ci = TRUE", {

  ts_ci <- get_drift_ts(res_tempssm, ci = TRUE)

  expect_s3_class(ts_ci, "ts")
  expect_equal(colnames(ts_ci), c("drift", "lwr", "upr"))
})


test_that("get_drift_ts checks inputs correctly", {

  expect_error(
    get_drift_ts(NULL),
    "`res` must be an object of class 'tempssm'"
  )

  expect_error(
    get_drift_ts(res_tempssm, ci = TRUE, ci_level = 1.2),
    "`ci_level` must be a numeric value between 0 and 1"
  )
})
