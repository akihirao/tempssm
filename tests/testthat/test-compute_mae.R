# test-compute_mae.R

test_that("compute_mae returns correct value", {

  y_true <- c(1, 2, 3)
  y_pred <- c(1, 3, 2)

  mae <- compute_mae(y_pred, y_true)

  expect_equal(mae, mean(c(0, 1, 1)))
})



test_that("compute_mae handles numeric values correctly", {

  y_true <- c(1.5, 2.5, 3.5)
  y_pred <- c(1.0, 3.0, 2.5)

  mae <- compute_mae(y_pred, y_true)

  expected <- mean(abs(y_true - y_pred))

  expect_equal(mae, expected)
})



test_that("compute_mae ignores NA values", {

  y_true <- c(1, NA, 3)
  y_pred <- c(1, 2, 2)

  mae <- compute_mae(y_pred, y_true)

  expected <- mean(abs(c(1 - 1, 3 - 2)))

  expect_equal(mae, expected)
})



test_that("compute_mae returns NA for NULL inputs", {

  expect_true(is.na(compute_mae(NULL, 1:3)))
  expect_true(is.na(compute_mae(1:3, NULL)))
})



test_that("length mismatch triggers error", {

  expect_error(
    compute_mae(1:3, 1:2),
    "must have the same length"
  )
})



test_that("compute_mae coerces inputs to numeric", {

  y_true <- c("1", "2", "3")
  y_pred <- c("1", "3", "2")

  mae <- compute_mae(y_pred, y_true)

  expect_equal(mae, mean(c(0, 1, 1)))
})



test_that("empty vectors return NA", {

  mae <- compute_mae(numeric(0), numeric(0))

  expect_true(is.na(mae) || is.nan(mae))
})

