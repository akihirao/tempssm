test_that(".prepare_tempssm_exogenous handles absent covariates", {
  expect_identical(
    .prepare_tempssm_exogenous(NULL),
    list(exo_names = NULL, exo_matrix = NULL)
  )
})


test_that(".prepare_tempssm_exogenous preserves covariate data and names", {
  exo_data <- ts(
    matrix(
      seq_len(16),
      ncol = 2,
      dimnames = list(NULL, c("pdo", "enso"))
    ),
    frequency = 4
  )
  withr::local_envvar(TEMPSSM_VERBOSITY = "")
  withr::local_options(list(tempssm.verbosity = "inform"))

  expect_silent(.prepare_tempssm_exogenous(exo_data))

  withr::local_options(list(tempssm.verbosity = "debug"))

  expect_message(
    exogenous <- .prepare_tempssm_exogenous(exo_data),
    "Including exogenous variables"
  )

  expect_identical(exogenous$exo_names, c("pdo", "enso"))
  expect_identical(exogenous$exo_matrix, as.matrix(exo_data))
  expect_identical(colnames(exogenous$exo_matrix), c("pdo", "enso"))
})
