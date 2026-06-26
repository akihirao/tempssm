test_that("model input checks reject complex time-series values", {
  complex_temp <- stats::ts(as.complex(1:4), frequency = 4)

  expect_error(
    .tempssm_prepare_model_inputs(complex_temp, na_action = "allow"),
    "temp_data.*numeric"
  )

  temp_ts <- stats::ts(1:4, frequency = 4)
  complex_exo <- stats::ts(as.complex(1:4), frequency = 4)
  complex_exo <- set_ts_name(complex_exo, label = "x", quiet = TRUE)

  expect_error(
    .tempssm_prepare_model_inputs(
      temp_data = temp_ts,
      exo_data = complex_exo
    ),
    "exo_data.*numeric"
  )
})


test_that("monthly conversion rejects unsupported Temp column types", {
  df_character <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = c("1", "2", "3")
  )
  df_complex <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = as.complex(1:3)
  )

  expect_error(
    convert_monthly_df_to_ts(df_character),
    "Temp.*numeric"
  )
  expect_error(
    convert_monthly_df_to_ts(df_complex),
    "Temp.*numeric"
  )
})


test_that("daily zoo conversion rejects unsupported variable types", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 3)
  zoo_character <- zoo::zoo(
    data.frame(Temp = c("1", "2", "3")),
    order.by = dates
  )
  zoo_complex <- zoo::zoo(
    data.frame(Temp = as.complex(1:3)),
    order.by = dates
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_character),
    "must be numeric"
  )
  expect_error(
    daily_zoo_to_monthly_ts(zoo_complex),
    "must be numeric"
  )
})


test_that("numeric scalar controls reject complex values", {
  temp_ts <- stats::ts(seq_len(12), frequency = 12)

  expect_error(
    tempssm(temp_ts, ar_order = 1 + 0i),
    "ar_order.*numeric"
  )

  expect_error(
    daily_zoo_to_monthly_ts(
      zoo::zoo(
        data.frame(Temp = seq_len(3)),
        order.by = as.Date("2001-01-01") + 0:2
      ),
      na_prop_max = 0.5 + 0i
    ),
    "na_prop_max.*numeric"
  )
})
