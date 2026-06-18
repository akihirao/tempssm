# test-convert_monthly_df_to_ts.R

test_that("convert_monthly_df_to_ts works for valid input", {
  df <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 12),
    Temp = rnorm(12)
  )

  ts_out <- convert_monthly_df_to_ts(df)

  expect_s3_class(ts_out, "ts")
  expect_identical(frequency(ts_out), 12)
  expect_identical(start(ts_out), c(2001, 1))
  expect_length(ts_out, 12)
})


test_that("warning is issued when months are missing", {
  df <- data.frame(
    Date = as.Date(c("2001-01-01", "2001-03-01")), # Feb missing
    Temp = c(10, 12)
  )

  expect_warning(
    convert_monthly_df_to_ts(df),
    "not strictly monthly"
  )
})


test_that("missing required columns triggers error", {
  df <- data.frame(
    Date = as.Date("2001-01-01")
  )

  expect_error(
    convert_monthly_df_to_ts(df),
    "must contain columns"
  )
})


test_that("non-Date column triggers error", {
  df <- data.frame(
    Date = "2001-01-01",
    Temp = 10
  )

  expect_error(
    convert_monthly_df_to_ts(df),
    "must be of class"
  )
})


test_that("non-data.frame input triggers error", {
  expect_error(
    convert_monthly_df_to_ts(NULL),
    "must be a data frame"
  )
})


test_that("start time is correctly set", {
  df <- data.frame(
    Date = seq(as.Date("1995-05-01"), by = "month", length.out = 6),
    Temp = rnorm(6)
  )

  ts_out <- convert_monthly_df_to_ts(df)

  expect_identical(start(ts_out), c(1995, 5))
})


test_that("NA values are allowed", {
  df <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = c(10, NA, 12)
  )

  ts_out <- convert_monthly_df_to_ts(df)

  expect_true(anyNA(ts_out))
})
