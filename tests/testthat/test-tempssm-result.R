test_that(".new_tempssm_result constructs fitted results", {
  temp_data <- ts(seq_len(8), frequency = 4)
  kfs <- list(
    alphahat = matrix(0, nrow = 8, ncol = 3),
    att = matrix(0, nrow = 8, ncol = 3)
  )
  state_names <- c("level", "slope", "arima1")

  res <- .new_tempssm_result(
    model = list(id = "model"),
    fit = list(id = "fit"),
    kfs = kfs,
    temp_data = temp_data,
    exogenous_data = NULL,
    ar_order = 1,
    use_season = FALSE,
    marginal = FALSE,
    model_call = quote(tempssm(temp_data = temp_data)),
    converged = TRUE,
    exo_names = NULL,
    state_names = state_names
  )

  expect_s3_class(res, "tempssm")
  expect_named(
    res,
    c(
      "model",
      "fit",
      "kfs",
      "temp_data",
      "exogenous_data",
      "ar_order",
      "use_season",
      "marginal",
      "call",
      "converged",
      "state_map"
    )
  )
  expect_identical(colnames(res$kfs$alphahat), state_names)
  expect_identical(colnames(res$kfs$att), state_names)
  expect_identical(res$state_map$all, state_names)
  expect_true(res$converged)
  expect_false(res$marginal)
})


test_that(".new_tempssm_result constructs fitting-error results", {
  temp_data <- ts(seq_len(8), frequency = 4)

  res <- .new_tempssm_result(
    model = NULL,
    fit = NULL,
    kfs = NULL,
    temp_data = temp_data,
    exogenous_data = NULL,
    ar_order = 1,
    use_season = TRUE,
    marginal = TRUE,
    model_call = quote(tempssm(temp_data = temp_data)),
    converged = FALSE,
    exo_names = NULL,
    state_names = character(0)
  )

  expect_s3_class(res, "tempssm")
  expect_null(res$model)
  expect_null(res$fit)
  expect_null(res$kfs)
  expect_false(res$converged)
  expect_true(res$marginal)
  expect_identical(res$temp_data, temp_data)
  expect_identical(res$state_map$all, character(0))
})
