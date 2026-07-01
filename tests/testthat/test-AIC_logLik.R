# test-AIC_logLik.R

test_that("AIC and logLik are consistent", {
  ll <- logLik(res_tempssm)
  expect_identical(
    get_aic(res_tempssm),
    -2 * as.numeric(ll) + 2 * attr(ll, "df")
  )
})


test_that("AIC.tempssm returns same value as get_aic", {
  a1 <- get_aic(res_tempssm)
  a2 <- AIC(res_tempssm)

  expect_identical(a1, a2)
})


test_that("AIC follows the selected likelihood setting", {
  ll <- logLik(res_tempssm, marginal = TRUE)
  expected <- -2 * as.numeric(ll) + 2 * attr(ll, "df")

  expect_identical(AIC(res_tempssm, marginal = TRUE), expected)
  expect_identical(get_aic(res_tempssm, marginal = TRUE), expected)
})
