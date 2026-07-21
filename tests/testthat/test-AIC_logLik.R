# test-AIC_logLik.R

test_that("AIC is intentionally unavailable for tempssm objects", {
  expect_error(
    AIC(res_tempssm),
    "AIC is not computed"
  )
})


test_that("get_aic is retained as a deprecated stopping wrapper", {
  expect_error(
    get_aic(res_tempssm),
    "deprecated"
  )
})
