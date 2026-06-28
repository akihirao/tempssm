test_that("monthly zoo conversion completes gaps and preserves metadata", {
  monthly_index <- zoo::as.yearmon(c(
    "2001-01", "2001-03", "2001-04"
  ))
  zoo_monthly <- zoo::zoo(
    matrix(c(1, 3, 4), ncol = 1),
    order.by = monthly_index
  )

  expect_warning(
    result <- .monthly_zoo_to_ts(zoo_monthly, var = "SST"),
    "Implicit missing months"
  )

  expect_s3_class(result, "ts")
  expect_identical(stats::start(result), c(2001, 1))
  expect_identical(stats::frequency(result), 12)
  expect_identical(colnames(result), "SST")
  expect_identical(as.numeric(result), c(1, NA, 3, 4))
})


test_that("monthly zoo conversion rejects undefined aggregated values", {
  zoo_monthly <- zoo::zoo(
    matrix(Inf, ncol = 1),
    order.by = zoo::as.yearmon("2001-01")
  )

  expect_error(
    .monthly_zoo_to_ts(zoo_monthly, var = "Temp"),
    "Inf"
  )
})
