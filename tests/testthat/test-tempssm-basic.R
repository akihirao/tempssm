# test-tempssm-basic.R

test_that("ssm runs and returns a tempssm object", {

  expect_s3_class(res_tempssm, "tempssm")
  expect_type(res_tempssm, "list")
})
