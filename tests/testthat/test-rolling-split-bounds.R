test_that(".make_rolling_split_bounds generates expanding bounds", {
  bounds <- .make_rolling_split_bounds(
    n_obs = 14,
    initial = 5,
    horizon = 3,
    step = 2,
    fixed_window = FALSE,
    allow_partial = FALSE
  )
  expected <- data.frame(
    train_start = c(1, 1, 1, 1),
    train_end = c(5, 7, 9, 11),
    test_start = c(6, 8, 10, 12),
    test_end = c(8, 10, 12, 14)
  )

  expect_identical(bounds, expected)
})


test_that(".make_rolling_split_bounds generates fixed-window bounds", {
  bounds <- .make_rolling_split_bounds(
    n_obs = 14,
    initial = 5,
    horizon = 3,
    step = 2,
    fixed_window = TRUE,
    allow_partial = FALSE
  )
  expected <- data.frame(
    train_start = c(1, 3, 5, 7),
    train_end = c(5, 7, 9, 11),
    test_start = c(6, 8, 10, 12),
    test_end = c(8, 10, 12, 14)
  )

  expect_identical(bounds, expected)
})


test_that(".make_rolling_split_bounds truncates a partial horizon", {
  complete <- .make_rolling_split_bounds(
    n_obs = 13,
    initial = 6,
    horizon = 4,
    step = 4,
    fixed_window = FALSE,
    allow_partial = FALSE
  )
  partial <- .make_rolling_split_bounds(
    n_obs = 13,
    initial = 6,
    horizon = 4,
    step = 4,
    fixed_window = FALSE,
    allow_partial = TRUE
  )

  expect_identical(complete$train_end, 6)
  expect_identical(partial$train_end, c(6, 10))
  expect_identical(partial$test_end, c(10, 13))
})


test_that(".make_rolling_split_bounds can return no complete bounds", {
  bounds <- .make_rolling_split_bounds(
    n_obs = 10,
    initial = 8,
    horizon = 4,
    step = 1,
    fixed_window = FALSE,
    allow_partial = FALSE
  )

  expect_s3_class(bounds, "data.frame")
  expect_named(
    bounds,
    c("train_start", "train_end", "test_start", "test_end")
  )
  expect_identical(nrow(bounds), 0L)
})
