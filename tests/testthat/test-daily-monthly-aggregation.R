test_that("daily month mean applies the missing-value threshold", {
  values <- c(1, 3, NA, 5)

  expect_identical(
    .daily_month_mean(values, na.rm = TRUE, na_prop_max = 0.25),
    3
  )
  expect_true(is.na(
    .daily_month_mean(values, na.rm = TRUE, na_prop_max = 0.24)
  ))
  expect_true(is.na(
    .daily_month_mean(values, na.rm = FALSE, na_prop_max = 1)
  ))
  expect_true(is.na(
    .daily_month_mean(rep(NA_real_, 4), na.rm = TRUE, na_prop_max = 1)
  ))
})


test_that("daily zoo aggregation returns calendar-month means", {
  dates <- as.Date(c(
    "2001-01-01", "2001-01-02",
    "2001-02-01", "2001-02-02"
  ))
  selected_zoo <- zoo::zoo(
    matrix(c(1, 3, 4, 8), ncol = 1),
    order.by = dates
  )

  result <- .aggregate_daily_zoo_monthly(
    selected_zoo,
    na.rm = TRUE,
    na_prop_max = 1
  )

  expect_s3_class(result, "zoo")
  expect_identical(NROW(result), 2L)
  expect_identical(as.numeric(zoo::coredata(result)), c(2, 6))
})
