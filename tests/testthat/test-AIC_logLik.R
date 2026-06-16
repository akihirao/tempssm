# test-AIC_logLik.R

test_that("AIC and logLik are consistent", {
  ll <- logLik(res_tempssm)
  expect_equal(
    get_aic(res_tempssm),
    -2 * as.numeric(ll) + 2 * attr(ll, "df")
  )
})


test_that("AIC.tempssm returns same value as get_aic", {
  a1 <- get_aic(res_tempssm)
  a2 <- AIC(res_tempssm)

  expect_equal(a1, a2)
})
