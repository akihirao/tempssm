test_that("compare_tempssm_aic returns a ranked comparison table", {
  models <- list(
    base = res_tempssm,
    exogenous = res_tempssm_exo
  )
  result <- compare_tempssm_aic(models)

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c(
      "model",
      "logLik",
      "df",
      "nobs",
      "observed_n",
      "start",
      "end",
      "frequency",
      "likelihood",
      "AIC",
      "delta_AIC",
      "weight"
    )
  )
  expected <- vapply(models, AIC, numeric(1))

  expect_identical(result$AIC, sort(result$AIC))
  expect_identical(result$AIC, unname(expected[result$model]))
  expect_identical(min(result$delta_AIC), 0)
  expect_lt(abs(sum(result$weight) - 1), sqrt(.Machine$double.eps))
  expect_true(all(result$likelihood == "marginal"))
  expect_true(all(result$nobs == length(temp_ts_test)))
  expect_true(all(result$observed_n == sum(!is.na(temp_ts_test))))
})


test_that("comparison validates list size and model names", {
  duplicated_names <- list(res_tempssm, res_tempssm_exo)
  names(duplicated_names) <- c("base", "base")

  expect_error(
    compare_tempssm_aic(list(base = res_tempssm)),
    "at least two models"
  )
  expect_error(
    compare_tempssm_aic(list(res_tempssm, res_tempssm_exo)),
    "non-empty name"
  )
  expect_error(
    compare_tempssm_aic(duplicated_names),
    "must be unique"
  )
})


test_that("comparison rejects invalid and non-converged models", {
  expect_error(
    compare_tempssm_aic(list(base = res_tempssm, invalid = list())),
    "Model.*invalid.*cannot be compared"
  )

  non_converged <- res_tempssm
  non_converged$converged <- FALSE
  expect_error(
    compare_tempssm_aic(
      list(base = res_tempssm, non_converged = non_converged)
    ),
    "Model.*non_converged.*did not converge"
  )
})


test_that("comparison requires identical response periods", {
  shifted <- res_tempssm
  shifted$temp_data <- ts(
    as.numeric(shifted$temp_data),
    start = c(2001, 1),
    frequency = 12
  )

  expect_error(
    compare_tempssm_aic(list(base = res_tempssm, shifted = shifted)),
    "same response period and time index"
  )
})


test_that("comparison requires identical response values and missingness", {
  changed <- res_tempssm
  changed$temp_data[1L] <- changed$temp_data[1L] + 1
  expect_error(
    compare_tempssm_aic(list(base = res_tempssm, changed = changed)),
    "same response values and missing-value pattern"
  )

  missing <- res_tempssm
  missing$temp_data[1L] <- NA_real_
  expect_error(
    compare_tempssm_aic(list(base = res_tempssm, missing = missing)),
    "same response values and missing-value pattern"
  )
})


test_that("comparison requires a common fitted likelihood type", {
  diffuse_model <- res_tempssm
  diffuse_model$marginal <- FALSE
  marginal_model <- res_tempssm_exo

  expect_error(
    compare_tempssm_aic(
      list(diffuse = diffuse_model, marginal = marginal_model)
    ),
    "diffuse and marginal likelihoods cannot be compared"
  )
})


test_that("comparison supports legacy diffuse-likelihood objects", {
  first <- res_tempssm
  second <- res_tempssm_exo
  first$marginal <- NULL
  second$marginal <- NULL

  result <- compare_tempssm_aic(list(first = first, second = second))

  expect_true(all(result$likelihood == "diffuse"))
})


test_that("compare_tempssm_aic does not modify fitted models", {
  models <- list(base = res_tempssm, exogenous = res_tempssm_exo)
  original <- models

  result <- compare_tempssm_aic(models)

  expect_s3_class(result, "tbl_df")
  expect_identical(models, original)
})
