#test-get_aic.R

test_that("get_aic returns a numeric value", {

  aic <- get_aic(res_tempssm)

  expect_type(aic, "double")
  expect_length(aic, 1)
})
