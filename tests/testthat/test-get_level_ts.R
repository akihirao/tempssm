#tests/testthat/test-get_level_ts.R

test_that("get_level_ts returns a ts object", {

  ts_obj <- get_level_ts(res_tempssm)

  expect_s3_class(ts_obj, "ts")
  expect_equal(NCOL(ts_obj), 1)
})


test_that("get_level_ts returns CI columns when ci = TRUE", {

  ts_ci <- get_level_ts(res_tempssm, ci = TRUE)

  expect_s3_class(ts_ci, "ts")
  expect_equal(colnames(ts_ci), c("ar1", "lwr", "upr"))
})



test_that("get_level_ts checks inputs correctly", {

  expect_error(
    get_level_ts(NULL),
    "`res` must be an object of class 'tempssm'"
  )
})
