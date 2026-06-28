test_that("daily zoo data preparation selects one numeric column", {
  dates <- as.Date("2001-01-01") + 0:3
  zoo_obj <- zoo::zoo(
    data.frame(Temp = 1:4, Salinity = 11:14),
    order.by = dates
  )

  prepared <- .prepare_daily_zoo_data(zoo_obj, var = "Salinity")

  expect_named(prepared, c("selected_zoo", "index", "values"))
  expect_s3_class(prepared$selected_zoo, "zoo")
  expect_identical(NCOL(prepared$selected_zoo), 1L)
  expect_identical(colnames(prepared$selected_zoo), "Salinity")
  expect_identical(prepared$index, dates)
  expect_identical(as.numeric(prepared$values), as.numeric(11:14))
})


test_that("daily zoo data preparation rejects invalid data", {
  dates <- as.Date("2001-01-01") + 0:2
  zoo_obj <- zoo::zoo(
    data.frame(Temp = letters[1:3]),
    order.by = dates
  )

  expect_error(
    .prepare_daily_zoo_data(1:3, var = "Temp"),
    "must be a <zoo> object"
  )
  expect_error(
    .prepare_daily_zoo_data(zoo_obj, var = "SST"),
    "not found"
  )
  expect_error(
    .prepare_daily_zoo_data(zoo_obj, var = "Temp"),
    "must be numeric"
  )
})
