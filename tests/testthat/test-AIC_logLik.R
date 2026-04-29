#test-AIC_logLik.R

test_that("AIC and logLik are consistent", {
  ll <- logLik(res_tempssm)
  expect_equal(
    get_aic(res_tempssm),
    -2 * as.numeric(ll) + 2 * attr(ll, "df")
  )
})
