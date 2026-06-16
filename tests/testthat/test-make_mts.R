# tests/testthat/test-make_mts.R

set.seed(123)

test_that(".make_mts converts ts to single-column mts", {
  y <- ts(1:10, start = c(2000, 1), frequency = 12)

  res <- .make_mts(y)

  expect_s3_class(res, "ts")
  expect_equal(NCOL(res), 1)
  expect_equal(colnames(res), "Temp")
})


test_that(".make_mts preserves time attributes", {
  y <- ts(1:10, start = c(2000, 5), frequency = 12)

  res <- .make_mts(y)

  expect_equal(start(res), start(y))
  expect_equal(frequency(res), frequency(y))
  expect_equal(time(res), time(y))
})


test_that(".make_mts preserves data values", {
  y <- ts(rnorm(10), frequency = 4)

  res <- .make_mts(y)

  expect_equal(as.numeric(res), as.numeric(y))
})


test_that(".make_mts works for matrix-like ts input", {
  y <- ts(matrix(1:10, ncol = 1), frequency = 12)

  res <- .make_mts(y)

  expect_equal(NCOL(res), 1)
  expect_equal(colnames(res), "Temp")
})


test_that(".make_mts handles NA values", {
  y <- ts(c(1, NA, 3, NA), frequency = 12)

  res <- .make_mts(y)

  expect_true(any(is.na(res)))
  expect_equal(as.numeric(res), as.numeric(y))
})
