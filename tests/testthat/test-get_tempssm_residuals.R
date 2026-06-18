# tests/testthat/test-get_tempssm_residuals.R

test_that("get_tempssm_residuals returns a numeric vector", {
  r <- get_tempssm_residuals(res_tempssm)

  expect_type(r, "double")
  expect_gt(length(r), 0)
  expect_true(all(is.finite(r)))
})


test_that("get_tempssm_residuals checks input class", {
  expect_error(
    get_tempssm_residuals("not a model"),
    "`res` must be an object of class <tempssm>."
  )
})
