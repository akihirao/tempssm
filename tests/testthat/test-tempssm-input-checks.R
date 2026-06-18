# tests/testthat/test-tempssm-input-checks.R

test_that(".tempssm_check_temp_ts accepts valid univariate ts", {
  temp_ts <- ts(rnorm(24), start = c(2000, 1), frequency = 12)

  expect_identical(.tempssm_check_temp_ts(temp_ts), temp_ts)
})


test_that(".tempssm_check_temp_ts rejects invalid temperature series", {
  expect_error(
    .tempssm_check_temp_ts(1:10),
    "must be a <ts> object"
  )

  expect_error(
    .tempssm_check_temp_ts(ts(rnorm(10), frequency = 1)),
    "frequency > 1"
  )

  multi_ts <- ts(matrix(rnorm(24), ncol = 2), frequency = 12)

  expect_error(
    .tempssm_check_temp_ts(multi_ts),
    "must be univariate"
  )
})


test_that(".tempssm_check_exo_ts accepts aligned named exogenous series", {
  temp_ts <- ts(rnorm(24), start = c(2000, 1), frequency = 12)
  exo_ts <- ts(
    matrix(rnorm(48), ncol = 2),
    start = c(2000, 1),
    frequency = 12
  )
  colnames(exo_ts) <- c("x1", "x2")

  expect_identical(.tempssm_check_exo_ts(temp_ts, exo_ts), exo_ts)
})


test_that(".tempssm_check_exo_ts rejects invalid exogenous series", {
  temp_ts <- ts(rnorm(24), start = c(2000, 1), frequency = 12)

  expect_error(
    .tempssm_check_exo_ts(temp_ts, 1:24),
    "must be a <ts> object"
  )

  exo_bad_freq <- ts(rnorm(8), start = c(2000, 1), frequency = 4)
  exo_bad_freq <- set_ts_name(exo_bad_freq, label = "x", quiet = TRUE)

  expect_error(
    .tempssm_check_exo_ts(temp_ts, exo_bad_freq),
    "Frequency"
  )

  exo_bad_time <- ts(rnorm(24), start = c(2001, 1), frequency = 12)
  exo_bad_time <- set_ts_name(exo_bad_time, label = "x", quiet = TRUE)

  expect_error(
    .tempssm_check_exo_ts(temp_ts, exo_bad_time),
    "Time index"
  )

  exo_no_name <- ts(rnorm(24), start = c(2000, 1), frequency = 12)

  expect_error(
    .tempssm_check_exo_ts(temp_ts, exo_no_name),
    "column name"
  )
})
