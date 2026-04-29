#tests/testthat/test-get_ar1_ts.R

test_that("get_ar1_ts returns a ts object", {

  ts_obj <- get_ar1_ts(res_tempssm)

  expect_s3_class(ts_obj, "ts")
  expect_equal(NCOL(ts_obj), 1)
})


test_that("get_ar1_ts returns CI columns when ci = TRUE", {

  ts_ci <- get_ar1_ts(res_tempssm, ci = TRUE)

  expect_s3_class(ts_ci, "ts")
  expect_equal(colnames(ts_ci), c("ar1", "lwr", "upr"))
})



test_that("get_ar1_ts checks inputs correctly", {

  expect_error(
    get_ar1_ts(NULL),
    "`res` must be an object of class 'tempssm'."
  )

  expect_error(
    get_ar1_ts(res_tempssm, ci = TRUE, ci_level = 1.2),
    "`ci_level` must be a numeric value between 0 and 1."
  )
})
