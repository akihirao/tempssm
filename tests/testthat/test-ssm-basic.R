# test-lgssm-basic.R

test_that("ssm runs and returns a tempssm object", {

  res <- ssm(temp_ts_test)

  expect_s3_class(res, "tempssm")
  expect_type(res, "list")
})
