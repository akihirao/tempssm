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


test_that("convert_monthly_df_to_ts accepts tibble input", {
  df <- tibble::tibble(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = c(10, 11, 12)
  )

  ts_out <- convert_monthly_df_to_ts(df)

  expect_s3_class(ts_out, "ts")
  expect_identical(start(ts_out), c(2001, 1))
  expect_identical(frequency(ts_out), 12)
  expect_identical(as.numeric(ts_out), c(10, 11, 12))
})


test_that("implicit missing months are inserted as explicit NA values", {
  df <- data.frame(
    Date = as.Date(c("2001-01-01", "2001-03-01")), # Feb missing
    Temp = c(10, 12)
  )

  expect_warning(
    ts_out <- convert_monthly_df_to_ts(df),
    "Implicit missing months"
  )

  expect_identical(start(ts_out), c(2001, 1))
  expect_length(ts_out, 3)
  expect_identical(as.numeric(ts_out), c(10, NA, 12))
})


test_that("unordered dates are caught and sorted before conversion", {
  df <- data.frame(
    Date = as.Date(c("2001-02-01", "2001-01-01")),
    Temp = c(12, 10)
  )

  expect_warning(
    ts_out <- convert_monthly_df_to_ts(df),
    "not ordered"
  )

  expect_identical(start(ts_out), c(2001, 1))
  expect_identical(as.numeric(ts_out), c(10, 12))
})


test_that("duplicate months trigger error", {
  df <- data.frame(
    Date = as.Date(c("2001-01-01", "2001-01-15")),
    Temp = c(10, 12)
  )

  expect_error(
    convert_monthly_df_to_ts(df),
    "duplicate months"
  )
})


test_that("missing dates trigger error", {
  df <- data.frame(
    Date = as.Date(c("2001-01-01", NA)),
    Temp = c(10, 12)
  )

  expect_error(
    convert_monthly_df_to_ts(df),
    "must not contain missing values"
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


test_that("non-numeric Temp column triggers error", {
  df <- data.frame(
    Date = as.Date("2001-01-01"),
    Temp = "10"
  )

  expect_error(
    convert_monthly_df_to_ts(df),
    "Temp.*numeric"
  )
})


test_that("base data frame list Temp column triggers informative error", {
  df <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 2),
    Temp = I(list(10, 11))
  )

  expect_error(
    convert_monthly_df_to_ts(df),
    "list column"
  )
})


test_that("tibble list Temp column triggers informative error", {
  df <- tibble::tibble(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 2),
    Temp = list(10, 11)
  )

  expect_error(
    convert_monthly_df_to_ts(df),
    "list column"
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


test_that("undefined Temp values trigger error", {
  df_nan <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = c(10, NaN, 12)
  )
  df_inf <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = c(10, Inf, 12)
  )

  expect_error(
    convert_monthly_df_to_ts(df_nan),
    "NaN"
  )

  expect_error(
    convert_monthly_df_to_ts(df_inf),
    "Inf"
  )
})


test_that("units Temp column is accepted and converted to numeric", {
  skip_if_not_installed("units")

  df <- data.frame(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = units::set_units(c(10, 11, 12), "K")
  )

  expect_warning(
    ts_out <- convert_monthly_df_to_ts(df),
    "converted to numeric"
  )

  expect_s3_class(ts_out, "ts")
  expect_false(inherits(ts_out, "units"))
  expect_identical(as.numeric(ts_out), c(10, 11, 12))
})


test_that("units Temp column in tibble input is converted to numeric", {
  skip_if_not_installed("units")

  df <- tibble::tibble(
    Date = seq(as.Date("2001-01-01"), by = "month", length.out = 3),
    Temp = units::set_units(c(10, 11, 12), "K")
  )

  expect_warning(
    ts_out <- convert_monthly_df_to_ts(df),
    "converted to numeric"
  )

  expect_s3_class(ts_out, "ts")
  expect_false(inherits(ts_out, "units"))
  expect_identical(as.numeric(ts_out), c(10, 11, 12))
})
