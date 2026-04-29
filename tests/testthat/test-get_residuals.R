#tests/testthat/test-get_residuals.R

test_that("get_residuals returns a numeric vector", {

  r <- get_residuals(res_tempssm)

  expect_type(r, "double")
  expect_true(length(r) > 0)
  expect_true(all(is.finite(r)))
})


test_that("get_residuals checks input class", {

  expect_error(
    get_residuals("not a model"),
    "`res` must be an object of class 'tempssm'."
  )
})
