test_that("autoplot.tempssm returns ggplot for each component", {
  components <- c("level", "drift", "season", "ar1")

  for (component in components) {
    p <- autoplot(
      res_tempssm,
      component = component,
      ci = FALSE
    )

    expect_s3_class(p, "ggplot")
  }
})


test_that("autoplot.tempssm passes extra arguments to component plotter", {
  p <- autoplot(
    res_tempssm,
    component = "level",
    ci = FALSE,
    ylab = "Sea surface temperature"
  )

  expect_identical(p$labels$y, "Sea surface temperature")
})


test_that("autoplot.tempssm returns gtable for all components", {
  g <- autoplot(res_tempssm, ci = FALSE)

  expect_s3_class(g, "gtable")
})


test_that("autoplot.tempssm propagates component input checks", {
  expect_error(
    autoplot(res_tempssm, component = "level", ci = "yes"),
    "`ci` must be a single logical value"
  )

  expect_error(
    autoplot(res_tempssm, component = "level", ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1."
  )
})

