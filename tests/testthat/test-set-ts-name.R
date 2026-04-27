test_that("set_ts_name works for univariate ts", {

  ts_uni <- ts(
    rnorm(120),
    start = c(2000, 1),
    frequency = 12
  )

  ts_named <- set_ts_name(ts_uni, label = "temperature")

  expect_s3_class(ts_named, "ts")
  expect_equal(NCOL(ts_named), 1)
  expect_equal(colnames(ts_named), "temperature")

  ## time attributes are preserved
  expect_equal(start(ts_named), start(ts_uni))
  expect_equal(frequency(ts_named), frequency(ts_uni))
  expect_equal(tsp(ts_named), tsp(ts_uni))
})


test_that("set_ts_name works for multivariate ts with matching labels", {

  ts_multi <- ts(
    matrix(rnorm(240), ncol = 3),
    start = c(1990, 1),
    frequency = 12
  )

  labels <- c("precip", "solar", "wind")

  ts_named <- set_ts_name(ts_multi, label = labels)

  expect_s3_class(ts_named, "ts")
  expect_equal(NCOL(ts_named), 3)
  expect_equal(colnames(ts_named), labels)

  ## time attributes are preserved
  expect_equal(tsp(ts_named), tsp(ts_multi))
})


test_that("single label is recycled for multivariate ts", {

  ts_multi <- ts(
    matrix(rnorm(120), ncol = 4),
    start = c(2010, 1),
    frequency = 12
  )

  ts_named <- set_ts_name(ts_multi, label = "exo")

  expect_equal(
    colnames(ts_named),
    rep("exo", 4)
  )
})


test_that("vector ts input is safely converted to single-column ts", {

  ts_vec <- ts(
    rnorm(60),
    start = c(2005, 1),
    frequency = 12
  )

  ts_named <- set_ts_name(ts_vec, label = "x")

  expect_s3_class(ts_named, "ts")
  expect_equal(NCOL(ts_named), 1)
  expect_equal(colnames(ts_named), "x")
})



test_that("non-ts input triggers error", {

  expect_error(
    set_ts_name(rnorm(10), label = "x"),
    "`ts_in` must be an object of class `ts`"
  )
})


test_that("non-character label triggers error", {

  ts_uni <- ts(rnorm(12), frequency = 12)

  expect_error(
    set_ts_name(ts_uni, label = 1),
    "`label` must be a character vector"
  )
})


test_that("invalid label length triggers error", {

  ts_multi <- ts(
    matrix(rnorm(100), ncol = 2),
    frequency = 12
  )

  expect_error(
    set_ts_name(ts_multi, label = c("a", "b", "c")),
    "Length of `label` must be 1 or equal"
  )
})





