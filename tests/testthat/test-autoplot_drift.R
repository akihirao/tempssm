test_that("autoplot_drift returns a ggplot object", {
  p1 <- autoplot_drift(res_tempssm)
  p2 <- autoplot_drift(res_tempssm, ci = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("autoplot_drift hides CI label in title by default", {
  p1 <- autoplot_drift(res_tempssm)
  p2 <- autoplot_drift(res_tempssm, show_ci_in_title = TRUE)

  expect_identical(p1$labels$title, "Drift component")
  expect_identical(p2$labels$title, "Drift component (95% CI)")
})


test_that("autoplot_drift uses compact rate label by default", {
  p <- autoplot_drift(res_tempssm)

  expect_identical(p$labels$y, "°C/yr")
})


test_that("autoplot_drift checks inputs correctly", {
  expect_error(
    autoplot_drift(NULL),
    "`res` must be an object of class <tempssm>."
  )

  expect_error(
    autoplot_drift(res_tempssm, ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1."
  )

  expect_error(
    autoplot_drift(res_tempssm, ci = "yes"),
    "`ci` must be a single logical value"
  )

  expect_error(
    autoplot_drift(res_tempssm, show_ci_in_title = "yes"),
    "`show_ci_in_title` must be a single logical value"
  )
})
