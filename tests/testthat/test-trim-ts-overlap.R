# tests/testthat/test-trim-ts-overlap.R

test_that("trim_ts_overlap works for valid input", {
  set.seed(123)

  exo_ts  <- ts(matrix(rnorm(120), ncol = 1),
                start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(
    temp_ts_test,
    exo_ts,
    temp_name = "temp",
    exo_name  = c("x1")
  )

  expect_true(is.list(res))
  expect_named(res, c("temperature", "exogenous"))

  expect_true(inherits(res$temperature, "ts"))
  expect_true(inherits(res$exogenous, "ts"))
})


test_that("trim_ts_overlap trims to overlapping period", {
  # length(temp_ts_test): 120
  exo_ts  <- ts(rnorm(50),  start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(temp_ts_test, exo_ts, exo_name = c("x1"))

  expect_true(length(res$temperature) == length(res$exogenous))
})


test_that("trim_ts_overlap assigns default exo names with warning", {
  exo_ts  <- ts(matrix(rnorm(120), ncol = 1),
                start = c(2000, 1), frequency = 12)

  expect_warning(
    res <- trim_ts_overlap(temp_ts_test, exo_ts)
  )

  expect_true(all(grepl("var", colnames(res$exogenous))))
})


test_that("trim_ts_overlap errors on incorrect exo_name length", {
  exo_ts  <- ts(matrix(rnorm(200), ncol = 2),
                start = c(2000, 1), frequency = 12)

  expect_error(
    trim_ts_overlap(temp_ts_test, exo_ts, exo_name = "x1")
  )
})



test_that("trim_ts_overlap handles multivariate exogenous series", {
  exo_ts  <- ts(matrix(rnorm(300), ncol = 3),
                start = c(2000, 1), frequency = 12)

  res <- trim_ts_overlap(
    temp_ts_test,
    exo_ts,
    exo_name = c("a", "b", "c")
  )

  expect_equal(colnames(res$exogenous), c("a", "b", "c"))
})



