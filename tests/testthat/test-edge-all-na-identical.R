test_that("model input checks reject all-NA temperature responses", {
  all_na_temp <- stats::ts(rep(NA_real_, 4), frequency = 4)

  expect_error(
    .tempssm_prepare_model_inputs(all_na_temp, na_action = "allow"),
    "at least one non-missing"
  )

  expect_error(
    .tempssm_prepare_model_inputs(all_na_temp, na_action = "warn"),
    "at least one non-missing"
  )
})


test_that("model input checks reject all-NA exogenous covariates", {
  temp_ts <- stats::ts(seq_len(4), frequency = 4)
  all_na_exo <- stats::ts(rep(NA_real_, 4), frequency = 4)
  all_na_exo <- set_ts_name(all_na_exo, label = "x", quiet = TRUE)

  expect_error(
    .tempssm_prepare_model_inputs(
      temp_data = temp_ts,
      exo_data = all_na_exo,
      na_action = "allow"
    ),
    "Exogenous covariates must be complete"
  )
})


test_that("conversion utilities preserve all-NA temperature fields", {
  monthly_df <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = rep(NA_real_, 3)
  )
  monthly_ts <- convert_monthly_df_to_ts(monthly_df)

  expect_true(all(is.na(monthly_ts)))

  daily_zoo <- zoo::zoo(
    data.frame(Temp = rep(NA_real_, 3)),
    order.by = as.Date("2001-01-01") + 0:2
  )

  expect_warning(
    daily_ts <- daily_zoo_to_monthly_ts(daily_zoo),
    "More than 30%"
  )
  expect_true(all(is.na(daily_ts)))
})


test_that("all-identical training series give undefined scaled errors", {
  constant_train <- stats::ts(rep(5, 10), frequency = 1)

  expect_identical(.scale_Q(constant_train, method = "naive"), 0)
  expect_true(
    is.na(.compute_mase(c(4, 6), c(5, 5), constant_train))
  )
})
