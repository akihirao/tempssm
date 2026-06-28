test_that("daily zoo controls preserve valid scalar values", {
  controls <- .prepare_daily_zoo_controls(
    var = "SST",
    na.rm = FALSE,
    na_prop_max = 0.25
  )

  expect_identical(
    controls,
    list(var = "SST", na.rm = FALSE, na_prop_max = 0.25)
  )
})


test_that("daily zoo controls reject undefined scalar values", {
  expect_error(
    .prepare_daily_zoo_controls(NA_character_, TRUE, 1),
    "var.*character scalar"
  )
  expect_error(
    .prepare_daily_zoo_controls("Temp", NA, 1),
    "na.rm.*logical scalar"
  )
  expect_error(
    .prepare_daily_zoo_controls("Temp", TRUE, Inf),
    "na_prop_max.*between 0 and 1"
  )
})
