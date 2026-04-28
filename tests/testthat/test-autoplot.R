test_that("autoplot.tempssm works for all components", {

  p1 <- autoplot(res_tempssm, component = "level")
  expect_s3_class(p1, "ggplot")

  p_all <- autoplot(res_tempssm)
  expect_true(inherits(p_all, "patchwork"))
})
