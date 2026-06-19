test_that("trim_prediction_intervals keeps intervals within width limit", {
  pred <- cbind(
    fit = c(1, 2, 3, 4),
    lwr = c(0, 0, 0, 0),
    upr = c(1, 2, 4, 5)
  )

  trimmed <- trim_prediction_intervals(pred, max_width = 2)

  expect_identical(nrow(trimmed), 2L)
  expect_identical(trimmed[, "fit"], pred[1:2, "fit"])
})


test_that("trim_prediction_intervals returns all rows when all widths pass", {
  pred <- cbind(
    fit = c(1, 2, 3),
    lwr = c(0, 1, 2),
    upr = c(1, 2, 3)
  )

  trimmed <- trim_prediction_intervals(pred, max_width = 1)

  expect_identical(trimmed, pred)
})


test_that("trim_prediction_intervals can return zero forecast rows", {
  pred <- cbind(
    fit = c(1, 2),
    lwr = c(0, 0),
    upr = c(2, 3)
  )

  trimmed <- trim_prediction_intervals(pred, max_width = 1)

  expect_identical(nrow(trimmed), 0L)
  expect_identical(colnames(trimmed), c("fit", "lwr", "upr"))
})


test_that("trim_prediction_intervals preserves prediction ts attributes", {
  pred <- ts(
    cbind(
      fit = c(1, 2, 3),
      lwr = c(0, 1, 1),
      upr = c(1, 2, 4)
    ),
    start = c(2020, 1),
    frequency = 12
  )

  trimmed <- trim_prediction_intervals(pred, max_width = 1)

  expect_true(inherits(trimmed, "ts"))
  expect_identical(stats::start(trimmed), c(2020, 1))
  expect_identical(stats::frequency(trimmed), 12)
  expect_identical(NROW(trimmed), 2L)
})


test_that("trim_prediction_intervals validates inputs", {
  pred <- cbind(fit = 1:2, lwr = 0:1, upr = 1:2)

  expect_error(
    trim_prediction_intervals(pred, max_width = NA_real_),
    "max_width"
  )

  expect_error(
    trim_prediction_intervals(1:3, max_width = 1),
    "pred"
  )

  expect_error(
    trim_prediction_intervals(pred[, "fit", drop = FALSE], max_width = 1),
    "lwr"
  )
})
