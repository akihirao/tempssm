# plot.Thermo

test_that("plot.tempssm runs without error", {

  res <- lgssm(temp_ts_test)

  expect_silent(plot(res))
  expect_silent(plot(res, ci = TRUE))

})
