# test-AIC_logLik.R

test_that("AIC is intentionally unavailable for tempssm objects", {
  expect_error(
    AIC(res_tempssm),
    "AIC is not computed"
  )
})
