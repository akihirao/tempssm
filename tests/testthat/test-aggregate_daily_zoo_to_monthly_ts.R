# tests/testthat/test-aggregate_daily_zoo_to_monthly_ts.R

test_that("aggregate_daily_zoo_to_monthly_ts works for valid input", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 60)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(60)),
    order.by = dates
  )

  res <- aggregate_daily_zoo_to_monthly_ts(zoo_obj)

  expect_s3_class(res, "ts")
  expect_identical(frequency(res), 12)
  expect_gt(length(res), 2) # 2 months
  expect_identical(colnames(res), "Temp")
})


test_that("var argument selects correct column", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(A = rnorm(30), B = rnorm(30)),
    order.by = dates
  )

  res <- aggregate_daily_zoo_to_monthly_ts(zoo_obj, var = "B")

  expect_identical(colnames(res), "B")
})


test_that("errors when variable not found", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(30)),
    order.by = dates
  )

  expect_error(
    aggregate_daily_zoo_to_monthly_ts(zoo_obj, var = "X"),
    "not found"
  )
})


test_that("errors for invalid na_prop_max", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(30)),
    order.by = dates
  )

  expect_error(aggregate_daily_zoo_to_monthly_ts(zoo_obj, na_prop_max = -1))
  expect_error(aggregate_daily_zoo_to_monthly_ts(zoo_obj, na_prop_max = 2))
})


test_that("errors when index is not Date/POSIXt", {
  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(10)),
    order.by = 1:10
  )

  expect_error(
    aggregate_daily_zoo_to_monthly_ts(zoo_obj),
    "Index of the zoo object"
  )
})


test_that("na.rm works correctly in aggregation", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  vals <- c(1, 2, NA, 3, 4, NA, 5, rep(1, 23))

  zoo_obj <- zoo::zoo(data.frame(Temp = vals), order.by = dates)

  res <- aggregate_daily_zoo_to_monthly_ts(zoo_obj, na.rm = TRUE)

  expect_false(is.na(res[1]))
})


test_that("warns when many NA in aggregated result", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 60)

  vals <- c(rep(NA, 40), rnorm(20)) # 66% NA

  zoo_obj <- zoo::zoo(data.frame(Temp = vals), order.by = dates)

  expect_warning(
    aggregate_daily_zoo_to_monthly_ts(zoo_obj),
    "More than 30%"
  )
})
