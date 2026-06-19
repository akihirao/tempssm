# tests/testthat/test-summary.R

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
  s <- summary(res_tempssm)
  ll <- logLik(res_tempssm)
  aic <- AIC(res_tempssm)

  expect_identical(s$logLik, as.numeric(ll))
  expect_identical(s$k, attr(ll, "df"))
  expect_identical(s$AIC, aic)
})


test_that("print.summary.tempssm works and returns object invisibly", {
  s <- summary(res_tempssm)

  expect_no_error(print(s))

  out <- print(s)
  expect_identical(out, s)
})


test_that("print.summary.tempssm outputs expected text", {
  s <- summary(res_tempssm)

  output <- capture.output(print(s))

  expect_true(any(grepl("tempssm summary", output, fixed = TRUE)))
  expect_true(any(grepl("Log-likelihood", output, fixed = TRUE)))
})


test_that("summary handles non-seasonal model", {
  data(niigata_sst)

  res <- tempssm(niigata_sst, use_season = FALSE, na_action = "allow")
  s <- summary(res)

  expect_true(is.na(s$variances$Q_season))
})


test_that("summary contains valid AR coefficients", {
  s <- summary(res_tempssm)

  expect_type(s$coef_ar, "list")
  expect_lte(s$coef_ar$AR_order, 1)
  expect_length(s$coef_ar$AR_coef, s$coef_ar$AR_order)
})


test_that("summary handles no exogenous variables", {
  s <- summary(res_tempssm)

  expect_true(is.null(s$exogenous) || is.character(s$exogenous))
})
