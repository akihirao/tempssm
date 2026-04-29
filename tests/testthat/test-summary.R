# test-summary.R

# test for summary object
test_that("summary.tempssm returns a valid summary object", {

  s <- summary(res_tempssm)

  expect_s3_class(s, "summary.tempssm")
  expect_type(s, "list")
})


# test for structure
test_that("summary contains expected components", {

  s <- summary(res_tempssm)

  expected_names <- c(
    "call",
    "logLik",
    "k",
    "AIC",
    "convergence",
    "variances",
    "coef_ar",
    "exogenous",
    "exogenous_coef"
  )

  expect_true(all(expected_names %in% names(s)))
})


test_that("summary uses logLik() and AIC() methods consistently", {

  s   <- summary(res_tempssm)
  ll  <- logLik(res_tempssm)
  aic <- AIC(res_tempssm)

  expect_equal(s$logLik, as.numeric(ll))
  expect_equal(s$k, attr(ll, "df"))
  expect_equal(s$AIC, aic)
})



