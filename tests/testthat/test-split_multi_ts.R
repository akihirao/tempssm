# tests/testthat/test-split_multi_ts.R

test_that("split_multi_ts splits multivariate ts correctly", {
  set.seed(123)

  multi_ts <- ts(
    matrix(rnorm(200), ncol = 2),
    start = c(2000, 1),
    frequency = 12
  )
  colnames(multi_ts) <- c("x1", "x2")

  res <- split_multi_ts(multi_ts)

  expect_true(is.list(res))
  expect_length(res, 2)

  expect_named(res, c("x1", "x2"))

  expect_true(all(sapply(res, function(x) inherits(x, "ts"))))
})


test_that("each split series is univariate ts", {
  multi_ts <- ts(
    matrix(rnorm(300), ncol = 3),
    start = c(2000, 1),
    frequency = 12
  )
  colnames(multi_ts) <- c("a", "b", "c")

  res <- split_multi_ts(multi_ts)

  expect_true(all(sapply(res, function(x) NCOL(x) == 1)))
})


test_that("split_multi_ts errors for non-ts input", {
  expect_error(split_multi_ts(NULL))
  expect_error(split_multi_ts("not ts"))
  expect_error(split_multi_ts(1:10))
})


test_that("split_multi_ts errors for univariate ts", {
  x <- ts(rnorm(100), start = c(2000, 1), frequency = 12)

  expect_error(split_multi_ts(x))
})


test_that("split_multi_ts preserves variable names", {
  multi_ts <- ts(
    matrix(rnorm(200), ncol = 2),
    start = c(2000, 1),
    frequency = 12
  )
  colnames(multi_ts) <- c("temp", "rain")

  res <- split_multi_ts(multi_ts)

  expect_equal(names(res), c("temp", "rain"))
})


test_that("each split series preserves time attributes", {
  multi_ts <- ts(
    matrix(rnorm(120), ncol = 2),
    start = c(1995, 1),
    frequency = 12
  )
  colnames(multi_ts) <- c("x1", "x2")

  res <- split_multi_ts(multi_ts)

  expect_equal(start(res[[1]]), start(multi_ts))
  expect_equal(frequency(res[[1]]), frequency(multi_ts))
})


test_that("each split ts has correct column name", {
  multi_ts <- ts(
    matrix(rnorm(120), ncol = 2),
    frequency = 12
  )
  colnames(multi_ts) <- c("var1", "var2")

  res <- split_multi_ts(multi_ts)

  expect_equal(colnames(res[[1]]), "var1")
  expect_equal(colnames(res[[2]]), "var2")
})


test_that("all elements are ts of length matching original", {
  multi_ts <- ts(
    matrix(rnorm(240), ncol = 3),
    frequency = 12
  )
  colnames(multi_ts) <- c("a", "b", "c")

  res <- split_multi_ts(multi_ts)

  expect_true(all(sapply(res, function(x) length(x) == nrow(multi_ts))))
})


test_that("order of variables is preserved", {
  multi_ts <- ts(
    matrix(rnorm(120), ncol = 3),
    frequency = 12
  )
  colnames(multi_ts) <- c("z", "x", "y")

  res <- split_multi_ts(multi_ts)

  expect_identical(names(res), c("z", "x", "y"))
})


test_that("handles more than 2 variables", {
  multi_ts <- ts(
    matrix(rnorm(300), ncol = 5),
    frequency = 12
  )
  colnames(multi_ts) <- paste0("v", 1:5)

  res <- split_multi_ts(multi_ts)

  expect_length(res, 5)
})
