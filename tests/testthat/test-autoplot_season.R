test_that("autoplot_season returns a ggplot object", {
  p1 <- autoplot_season(res_tempssm)
  p2 <- autoplot_season(res_tempssm, ci = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("autoplot_season hides CI label in title by default", {
  p1 <- autoplot_season(res_tempssm)
  p2 <- autoplot_season(res_tempssm, show_ci_in_title = TRUE)

  expect_identical(p1$labels$title, "Seasonal component")
  expect_identical(p2$labels$title, "Seasonal component (95% CI)")
})


test_that("autoplot_season checks inputs correctly", {
  expect_error(
    autoplot_season(NULL),
    "`res` must be an object of class <tempssm>."
  )

  expect_error(
    autoplot_season(res_tempssm, ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1."
  )

  expect_error(
    autoplot_season(res_tempssm, ci = "yes"),
    "`ci` must be a single logical value"
  )

  expect_error(
    autoplot_season(res_tempssm, show_ci_in_title = "yes"),
    "`show_ci_in_title` must be a single logical value"
  )
})
