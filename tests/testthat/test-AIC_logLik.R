# test-AIC_logLik.R

test_that("AIC methods are consistent with the selected log-likelihood", {
  ll <- logLik(res_tempssm)
  expected <- -2 * as.numeric(ll) + 2 * attr(ll, "df")

  expect_identical(get_aic(res_tempssm), expected)
  expect_identical(AIC(res_tempssm), expected)
})


test_that("AIC methods support an explicit diffuse likelihood", {
  ll <- logLik(res_tempssm, marginal = FALSE)
  expected <- -2 * as.numeric(ll) + 2 * attr(ll, "df")

  expect_identical(AIC(res_tempssm, marginal = FALSE), expected)
  expect_identical(get_aic(res_tempssm, marginal = FALSE), expected)
})
