test_that("get_aic is retained as deprecated unavailable API", {
  expect_error(
    get_aic(res_tempssm),
    "deprecated"
  )
})
