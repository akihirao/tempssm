# tests/testthat/test-handle_non_convergence.R

test_that(".handle_non_convergence handles NULL result", {
  y <- ts(1:10, frequency = 12)
  fold <- list(fold = 1)

  expect_warning(
    res <- .handle_non_convergence(NULL, fold, y, y),
    regexp = "did not converge"
  )

  expect_type(res, "list")
  expect_false(res$converged)
  expect_null(res$y_pred)
  expect_null(res$model)
})


test_that(".handle_non_convergence handles non-converged model", {
  y <- ts(1:10, frequency = 12)
  fold <- list(fold = 2)

  mock_res <- list(
    converged = FALSE,
    model = "dummy_model"
  )

  expect_warning(
    res <- .handle_non_convergence(mock_res, fold, y, y),
    regexp = "did not converge"
  )

  expect_false(res$converged)
  expect_equal(res$model, "dummy_model")
})


test_that(".handle_non_convergence returns NULL when converged", {
  y <- ts(1:10, frequency = 12)
  fold <- list(fold = 3)

  mock_res <- list(
    converged = TRUE,
    model = "ok"
  )

  res <- .handle_non_convergence(mock_res, fold, y, y)

  expect_null(res)
})


test_that("warning includes fold id", {
  y <- ts(1:10, frequency = 12)
  fold <- list(fold = 999)

  expect_warning(
    .handle_non_convergence(NULL, fold, y, y),
    regexp = "999"
  )
})


test_that("returned structure is consistent", {
  y <- ts(1:5, frequency = 4)
  fold <- list(fold = 10)

  res <- suppressWarnings(
    .handle_non_convergence(NULL, fold, y, y)
  )

  expect_named(res, c(
    "fold", "converged", "y_train", "y_test", "y_pred", "model"
  ))

  expect_identical(res$y_train, y)
  expect_identical(res$y_test, y)
})
