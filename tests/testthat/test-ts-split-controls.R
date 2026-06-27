test_that(".prepare_ts_split_controls preserves valid controls", {
  controls <- .prepare_ts_split_controls(
    initial = 24,
    horizon = 6,
    step = 12,
    fixed_window = TRUE,
    allow_partial = FALSE,
    n_obs = 48
  )

  expect_identical(
    controls,
    list(
      initial = 24,
      horizon = 6,
      step = 12,
      fixed_window = TRUE,
      allow_partial = FALSE
    )
  )
})


test_that(".prepare_ts_split_controls validates count lengths and types", {
  expect_error(
    .prepare_ts_split_controls(
      c(12, 24), 6, 12, FALSE, FALSE, 48
    ),
    "initial.*length one"
  )
  expect_error(
    .prepare_ts_split_controls(
      24, "6", 12, FALSE, FALSE, 48
    ),
    "horizon.*numeric"
  )
})


test_that(".prepare_ts_split_controls rejects invalid count values", {
  expect_error(
    .prepare_ts_split_controls(0, 6, 12, FALSE, FALSE, 48),
    "must be positive integers"
  )
  expect_error(
    .prepare_ts_split_controls(24, Inf, 12, FALSE, FALSE, 48),
    "must be positive integers"
  )
  expect_error(
    .prepare_ts_split_controls(24, 6, 1.5, FALSE, FALSE, 48),
    "must be positive integers"
  )
  expect_error(
    .prepare_ts_split_controls(NA_real_, 6, 12, FALSE, FALSE, 48),
    "must be positive integers"
  )
})


test_that(".prepare_ts_split_controls validates logical controls", {
  expect_error(
    .prepare_ts_split_controls(24, 6, 12, NA, FALSE, 48),
    "fixed_window.*logical scalar"
  )
  expect_error(
    .prepare_ts_split_controls(24, 6, 12, FALSE, NA, 48),
    "allow_partial.*logical scalar"
  )
  expect_error(
    .prepare_ts_split_controls(24, 6, 12, 1, FALSE, 48),
    "fixed_window.*logical"
  )
})


test_that(".prepare_ts_split_controls requires room after initial window", {
  expect_error(
    .prepare_ts_split_controls(48, 6, 12, FALSE, FALSE, 48),
    "initial.*smaller than the length"
  )
})
