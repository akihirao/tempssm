test_that("autoplot_drift returns a ggplot object", {

  p1 <- autoplot_drift(res_tempssm)
  p2 <- autoplot_drift(res_tempssm, ci = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("autoplot_drift checks inputs correctly", {

  expect_error(
    autoplot_drift(NULL),
    "`res` must be an object returned by tempssm"
  )

  expect_error(
    autoplot_drift(res_tempssm, ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1"
  )

  expect_error(
    autoplot_drift(res_tempssm, ci = "yes"),
    "`ci` must be a single logical value"
  )
})

