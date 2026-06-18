# tests/testthat/test-plot_tempssm_residual_diagnostics.R

test_that("plot_tempssm_residual_diagnostics runs without error", {
  expect_silent(
    plot_tempssm_residual_diagnostics(res_tempssm, save = FALSE)
  )
})
