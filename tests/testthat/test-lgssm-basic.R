# test-lgssm-basic.R

test_that("lgssm runs and returns a ThermoSSM object", {

  res <- lgssm(temp_ts_test)

  expect_s3_class(res, "ThermoSSM")
  expect_type(res, "list")
})
