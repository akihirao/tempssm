test_that("autoplot_level returns a ggplot object", {
  p1 <- autoplot_level(res_tempssm)
  p2 <- autoplot_level(res_tempssm, ci = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("autoplot_level adds confidence interval layers when requested", {
  with_ci <- autoplot_level(res_tempssm, ci = TRUE)
  without_ci <- autoplot_level(res_tempssm, ci = FALSE)

  expect_length(with_ci$layers, 2L)
  expect_s3_class(with_ci$layers[[1]]$geom, "GeomLine")
  expect_s3_class(with_ci$layers[[2]]$geom, "GeomRibbon")
  expect_length(without_ci$layers, 1L)
  expect_s3_class(without_ci$layers[[1]]$geom, "GeomLine")
})


test_that("autoplot_level hides CI label in title by default", {
  p1 <- autoplot_level(res_tempssm)
  p2 <- autoplot_level(res_tempssm, show_ci_in_title = TRUE)

  expect_identical(p1$labels$title, "Level component")
  expect_identical(p2$labels$title, "Level component (95% CI)")
})


test_that("autoplot_level uses compact unit label by default", {
  p <- autoplot_level(res_tempssm)

  expect_identical(rlang::as_label(p$mapping$x), "time")
  expect_identical(p$labels$x, "Time (year)")
  expect_identical(p$labels$y, expression(Temp. ~ (degree * C)))
})


test_that("autoplot_level checks inputs correctly", {
  expect_error(
    autoplot_level(NULL),
    "`res` must be an object of class <tempssm>."
  )

  expect_error(
    autoplot_level(res_tempssm, ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1."
  )

  expect_error(
    autoplot_level(res_tempssm, ci = "yes"),
    "`ci` must be a single logical value"
  )

  expect_error(
    autoplot_level(res_tempssm, show_ci_in_title = "yes"),
    "`show_ci_in_title` must be a single logical value"
  )
})
