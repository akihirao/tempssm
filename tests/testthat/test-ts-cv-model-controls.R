test_that("CV model controls preserve valid values", {
  controls <- .prepare_ts_cv_model_controls(
    ar_order = 2,
    use_season = FALSE
  )

  expect_named(controls, c("ar_order", "use_season"))
  expect_identical(controls$ar_order, 2)
  expect_false(controls$use_season)
})


test_that("CV model controls do not issue model-fitting warnings", {
  expect_no_warning(
    controls <- .prepare_ts_cv_model_controls(
      ar_order = 5,
      use_season = TRUE
    )
  )

  expect_identical(controls$ar_order, 5)
})


test_that("CV model controls reject invalid autoregressive orders", {
  expect_error(
    .prepare_ts_cv_model_controls(0, TRUE),
    "ar_order.*integer >= 1"
  )
  expect_error(
    .prepare_ts_cv_model_controls(1.5, TRUE),
    "ar_order.*integer >= 1"
  )
  expect_error(
    .prepare_ts_cv_model_controls(NA_real_, TRUE),
    "ar_order.*integer >= 1"
  )
  expect_error(
    .prepare_ts_cv_model_controls(c(1, 2), TRUE),
    "ar_order.*length one"
  )
})


test_that("CV model controls reject invalid seasonal controls", {
  expect_error(
    .prepare_ts_cv_model_controls(1, NA),
    "use_season.*logical scalar"
  )
  expect_error(
    .prepare_ts_cv_model_controls(1, c(TRUE, FALSE)),
    "use_season.*length one"
  )
})
