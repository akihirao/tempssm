#tests/testthat/test-plot_residual_diagnostics.R

test_that("plot_residual_diagnostics runs without error", {

  expect_silent(
    plot_residual_diagnostics(res_tempssm, save = FALSE)
  )
})
