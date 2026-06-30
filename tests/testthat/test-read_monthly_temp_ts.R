# tests/testthat/test-read_monthly_temp_ts.R

test_that("read_monthly_temp_ts works for valid CSV", {
  tmp <- tempfile(fileext = ".csv")

  writeLines(
    c(
      "Year,Month,Temp",
      "2001,1,10.4",
      "2001,2,8.2",
      "2001,3,NA",
      "2001,4,13.6"
    ),
    tmp
  )

  ts_out <- read_monthly_temp_ts(tmp)

  expect_s3_class(ts_out, "ts")
  expect_identical(frequency(ts_out), 12)
  expect_identical(start(ts_out), c(2001, 1))
  expect_length(ts_out, 4)
})


test_that("errors when csv is not character scalar", {
  expect_error(
    read_monthly_temp_ts(123),
    "`csv` must be a file path"
  )
  expect_error(
    read_monthly_temp_ts(c("a", "b")),
    "`csv` must be a file path"
  )
})


test_that("errors when file does not exist", {
  expect_error(
    read_monthly_temp_ts("nonexistent_file.csv"),
    "File does not exist"
  )
})


test_that("errors when required columns are missing", {
  tmp <- tempfile(fileext = ".csv")

  writeLines(
    c(
      "Year,Month",
      "2001,1",
      "2001,2"
    ),
    tmp
  )

  expect_error(
    read_monthly_temp_ts(tmp),
    "must contain columns"
  )
})


test_that("errors when CSV has no rows", {
  tmp <- tempfile(fileext = ".csv")

  writeLines(
    "Year,Month,Temp",
    tmp
  )

  expect_error(
    read_monthly_temp_ts(tmp),
    "contains no rows"
  )
})


test_that("errors when month values are invalid", {
  tmp <- tempfile(fileext = ".csv")

  writeLines(
    c(
      "Year,Month,Temp",
      "2001,13,10.4",
      "2001,2,8.2"
    ),
    tmp
  )

  expect_error(
    read_monthly_temp_ts(tmp),
    "invalid month"
  )
})


test_that("errors when year-month rows are duplicated", {
  tmp <- tempfile(fileext = ".csv")

  writeLines(
    c(
      "Year,Month,Temp",
      "2001,1,10.4",
      "2001,1,8.2"
    ),
    tmp
  )

  expect_error(
    read_monthly_temp_ts(tmp),
    "duplicate year-month"
  )
})


test_that("errors when year-month rows are not increasing", {
  tmp <- tempfile(fileext = ".csv")

  writeLines(
    c(
      "Year,Month,Temp",
      "2001,2,8.2",
      "2001,1,10.4"
    ),
    tmp
  )

  expect_error(
    read_monthly_temp_ts(tmp),
    "strictly increasing"
  )
})


test_that("implicit missing months in CSV input are explicit NA values", {
  tmp <- tempfile(fileext = ".csv")

  writeLines(
    c(
      "Year,Month,Temp",
      "2001,1,10.4",
      "2001,3,13.6"
    ),
    tmp
  )

  expect_warning(
    ts_out <- read_monthly_temp_ts(tmp),
    "Implicit missing months"
  )

  expect_identical(start(ts_out), c(2001, 1))
  expect_length(ts_out, 3)
  expect_identical(as.numeric(ts_out), c(10.4, NA, 13.6))
})


test_that("missing Year or Month values trigger error", {
  rows <- c(
    "NA,1,10.4",
    "2001,NA,10.4"
  )

  for (row in rows) {
    tmp <- tempfile(fileext = ".csv")
    writeLines(c("Year,Month,Temp", row), tmp)

    expect_error(
      read_monthly_temp_ts(tmp),
      "must not contain missing values"
    )
  }
})


test_that("Year and Month columns must contain whole numeric values", {
  cases <- list(
    list(row = "year,1,10.4", message = "must be numeric"),
    list(row = "2001,month,10.4", message = "must be numeric"),
    list(row = "2001.5,1,10.4", message = "whole numbers"),
    list(row = "2001,1.5,10.4", message = "whole numbers")
  )

  for (case in cases) {
    tmp <- tempfile(fileext = ".csv")
    writeLines(c("Year,Month,Temp", case$row), tmp)

    expect_error(
      read_monthly_temp_ts(tmp),
      case$message
    )
  }
})


test_that("Temp column must be numeric", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("Year,Month,Temp", "2001,1,warm"), tmp)

  expect_error(
    read_monthly_temp_ts(tmp),
    "Temp.*numeric"
  )
})


test_that("undefined Temp values trigger error", {
  cases <- c("NaN", "Inf", "-Inf")

  for (value in cases) {
    tmp <- tempfile(fileext = ".csv")
    writeLines(c("Year,Month,Temp", paste("2001,1", value, sep = ",")), tmp)

    expect_error(
      read_monthly_temp_ts(tmp),
      "NaN|Inf"
    )
  }
})
