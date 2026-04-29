#test-AIC.R

test_that("AIC.tempssm returns same value as get_aic", {

  a1 <- get_aic(res_tempssm)
  a2 <- AIC(res_tempssm)

  expect_equal(a1, a2)
})
