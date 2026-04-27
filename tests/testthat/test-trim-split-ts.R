test_that("trim_ts_overlap and split_multi_ts work together", {

  ## setup
  temp_ts <- ts(
    rnorm(120),
    start = c(2000, 1),
    frequency = 12
  )

  exo_ts <- ts(
    matrix(rnorm(240), ncol = 2),
    start = c(2001, 1),
    frequency = 12
  )
  colnames(exo_ts) <- c("precip", "solar")

  ## run
  trimmed <- trim_ts_overlap(
    temp_ts,
    exo_ts,
    temp_name = "temp",
    exo_name  = c("precip", "solar")
  )

  split_exo <- split_multi_ts(trimmed$exogenous)

  ## expectations
  expect_type(trimmed, "list")
  expect_named(trimmed, c("temperature", "exogenous"))

  expect_s3_class(trimmed$temperature, "ts")
  expect_s3_class(trimmed$exogenous, "ts")

  expect_equal(colnames(trimmed$temperature), "temp")
  expect_equal(colnames(trimmed$exogenous), c("precip", "solar"))

  expect_length(split_exo, 2)
  expect_named(split_exo, c("precip", "solar"))

  expect_s3_class(split_exo[[1]], "ts")
  expect_s3_class(split_exo[[2]], "ts")

  expect_equal(colnames(split_exo[[1]]), "precip")
  expect_equal(colnames(split_exo[[2]]), "solar")

})

