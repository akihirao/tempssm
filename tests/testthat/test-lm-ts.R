# test-lm-ts.R

test_that("lm_ts returns lm object for valid input", {

  fit <- lm_ts(temp_ts_test)

  expect_s3_class(fit, "lm")
})


test_that("lm_ts produces valid coefficients", {

  fit <- lm_ts(temp_ts_test)

  coef <- coefficients(fit)

  expect_true(length(coef) == 2)   # intercept + slope
  expect_true(all(is.finite(coef)))
})


test_that("lm_ts errors on non-ts input", {

  expect_error(lm_ts(NULL))
  expect_error(lm_ts(1:10))
  expect_error(lm_ts("not ts"))
})


test_that("lm_ts errors on multivariate ts input", {

  multi_ts <- ts(
    matrix(rnorm(200), ncol = 2),
    start = c(2000, 1),
    frequency = 12
  )

  expect_error(lm_ts(multi_ts))
})

