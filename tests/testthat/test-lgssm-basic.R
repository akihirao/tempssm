# test-lgssm-basic.R

test_that("lgssm runs and returns a tempssm object", {

  res <- lgssm(temp_ts_test)

  expect_s3_class(res, "tempssm")
  expect_type(res, "list")
})
