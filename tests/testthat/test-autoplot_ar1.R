test_that("autoplot_ar1 returns a ggplot object", {

  p1 <- autoplot_ar1(res_tempssm)
  p2 <- autoplot_ar1(res_tempssm, ci = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("autoplot_ar1 checks inputs correctly", {

  expect_error(
    autoplot_ar1(NULL),
    "`res` must be an object returned by tempssm"
  )

  expect_error(
    autoplot_ar1(res_tempssm, ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1"
  )

  expect_error(
    autoplot_ar1(res_tempssm, ci = "yes"),
    "`ci` must be a single logical value"
  )
})
