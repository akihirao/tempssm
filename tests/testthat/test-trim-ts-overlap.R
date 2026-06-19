# tests/testthat/test-trim-ts-overlap.R

test_that("trim_ts_overlap works for valid input", {
  set.seed(123)

  exo_ts <- ts(matrix(rnorm(120), ncol = 1),
    start = c(2000, 1), frequency = 12
  )

  res <- trim_ts_overlap(
    temp_ts_test,
    exo_ts,
    temp_name = "temp",
    exo_name  = "x1"
  )

  expect_type(res, "list")
  expect_named(res, c("temperature", "exogenous"))

  expect_s3_class(res$temperature, "ts")
  expect_s3_class(res$exogenous, "ts")
})


test_that("trim_ts_overlap accepts explicit temp_data and exo_data", {
  exo_ts <- ts(rnorm(120), start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(
    temp_data = temp_ts_test,
    exo_data = exo_ts,
    exo_name = "x1"
  )

  expect_s3_class(res$temperature, "ts")
  expect_s3_class(res$exogenous, "ts")
})


test_that("trim_ts_overlap accepts legacy temp_ts and exo_ts aliases", {
  exo_ts <- ts(rnorm(120), start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(
    temp_ts = temp_ts_test,
    exo_ts = exo_ts,
    exo_name = "x1"
  )

  expect_s3_class(res$temperature, "ts")
  expect_s3_class(res$exogenous, "ts")
})


test_that("trim_ts_overlap rejects mixed current and legacy inputs", {
  exo_ts <- ts(rnorm(120), start = c(2000, 1), frequency = 12)

  expect_error(
    trim_ts_overlap(
      temp_data = temp_ts_test,
      exo_data = exo_ts,
      temp_ts = temp_ts_test,
      exo_name = "x1"
    ),
    "Use either"
  )

  expect_error(
    trim_ts_overlap(
      temp_data = temp_ts_test,
      exo_data = exo_ts,
      exo_ts = exo_ts,
      exo_name = "x1"
    ),
    "Use either"
  )
})


test_that("trim_ts_overlap trims to overlapping period", {
  # length(temp_ts_test): 120
  exo_ts <- ts(rnorm(50), start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(temp_ts_test, exo_ts, exo_name = "x1")

  expect_length(res$temperature, length(res$exogenous))
})


test_that("trim_ts_overlap assigns default exo names with warning", {
  exo_ts <- ts(matrix(rnorm(120), ncol = 1),
    start = c(2000, 1), frequency = 12
  )

  expect_warning(
    res <- trim_ts_overlap(temp_ts_test, exo_ts)
  )

  expect_true(all(grepl("var", colnames(res$exogenous), fixed = TRUE)))
})


test_that("trim_ts_overlap errors on incorrect exo_name length", {
  exo_ts <- ts(matrix(rnorm(200), ncol = 2),
    start = c(2000, 1), frequency = 12
  )

  expect_error(
    trim_ts_overlap(temp_ts_test, exo_ts, exo_name = "x1")
  )
})


test_that("trim_ts_overlap handles multivariate exogenous series", {
  exo_ts <- ts(matrix(rnorm(300), ncol = 3),
    start = c(2000, 1), frequency = 12
  )

  res <- trim_ts_overlap(
    temp_ts_test,
    exo_ts,
    exo_name = c("a", "b", "c")
  )

  expect_identical(colnames(res$exogenous), c("a", "b", "c"))
})


test_that("errors when temp_ts is not ts", {
  exo_ts <- ts(rnorm(12), frequency = 12)

  expect_error(
    trim_ts_overlap(rnorm(12), exo_ts),
    "must be an object of class"
  )
})

test_that("errors when exo_ts is not ts", {
  expect_error(
    trim_ts_overlap(temp_ts_test, rnorm(12)),
    "must be an object of class"
  )
})


test_that("temp_name is correctly assigned", {
  exo_ts <- ts(rnorm(120), start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(
    temp_ts_test,
    exo_ts,
    temp_name = "temperature_test",
    exo_name = "x"
  )

  expect_identical(colnames(res$temperature), "temperature_test")
})


test_that("overlap period is correctly computed", {
  temp_ts <- ts(rnorm(120), start = c(2000, 1), frequency = 12)
  exo_ts <- ts(rnorm(60), start = c(2005, 1), frequency = 12)

  res <- trim_ts_overlap(temp_ts, exo_ts, exo_name = "x")

  expect_identical(start(res$temperature), c(2005, 1))
  expect_identical(end(res$temperature), end(exo_ts))
})


test_that("single exogenous variable uses set_ts_name", {
  exo_ts <- ts(rnorm(120), start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(
    temp_ts_test,
    exo_ts,
    exo_name = "x"
  )

  expect_identical(colnames(res$exogenous), "x")
})


test_that("multivariate exogenous uses direct colnames assignment", {
  exo_ts <- ts(matrix(rnorm(240), ncol = 2),
    start = c(2000, 1), frequency = 12
  )

  res <- trim_ts_overlap(
    temp_ts_test,
    exo_ts,
    exo_name = c("a", "b")
  )

  expect_identical(colnames(res$exogenous), c("a", "b"))
})
