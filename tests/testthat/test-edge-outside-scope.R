test_that("model input checks reject multivariate temperature responses", {
  temp_multi <- stats::ts(
    matrix(seq_len(24), ncol = 2),
    frequency = 4
  )

  expect_error(
    .tempssm_prepare_model_inputs(temp_multi),
    "temp_data.*univariate"
  )
})


test_that("model input checks reject more exogenous fields than observations", {
  temp_ts <- stats::ts(seq_len(4), frequency = 4)
  exo_ts <- stats::ts(
    matrix(seq_len(20), nrow = 4, ncol = 5),
    frequency = 4
  )
  colnames(exo_ts) <- paste0("x", seq_len(NCOL(exo_ts)))

  expect_error(
    .tempssm_prepare_model_inputs(
      temp_data = temp_ts,
      exo_data = exo_ts
    ),
    "exogenous variables.*observations"
  )
})
