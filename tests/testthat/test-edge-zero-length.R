test_that("monthly data-frame conversion rejects zero-length input", {
  empty_df <- data.frame(
    Date = as.Date(character(0)),
    Temp = numeric(0)
  )

  expect_error(
    convert_monthly_df_to_ts(empty_df),
    "No observations are available"
  )
})


test_that("daily zoo conversion rejects zero-length input", {
  empty_zoo <- zoo::zoo(
    data.frame(Temp = numeric(0)),
    order.by = as.Date(character(0))
  )

  expect_error(
    daily_zoo_to_monthly_ts(empty_zoo),
    "No observations are available"
  )
})


test_that("MAE handles zero-length vectors explicitly", {
  expect_true(is.na(.compute_mae(numeric(0), numeric(0))))
})
