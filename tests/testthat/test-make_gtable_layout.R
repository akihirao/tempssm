# test-make_gtable_layout

test_that("grid_arrange_base works with multiple ggplots", {
  skip_if_not_installed("ggplot2")

  p1 <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point()

  p2 <- ggplot2::ggplot(mtcars, ggplot2::aes(factor(cyl), mpg)) +
    ggplot2::geom_boxplot()

  plots <- list(p1, p2)

  expect_silent(
    .make_gtable_layout(plots)
  )
})

test_that("grid_arrange_base respects nrow/ncol", {
  skip_if_not_installed("ggplot2")

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point()

  plots <- list(p, p, p, p)

  expect_silent(
    .make_gtable_layout(plots, nrow = 2, ncol = 2)
  )
})

test_that("grid_arrange_base handles single plot", {
  skip_if_not_installed("ggplot2")

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point()

  expect_silent(
    .make_gtable_layout(list(p))
  )
})

test_that(".make_gtable_layout requires a list of plots", {
  expect_error(
    .make_gtable_layout(NULL),
    "`plots` must be a list."
  )
})

test_that("autoplot.tempssm validates component", {
  expect_error(
    autoplot(res_tempssm, component = c("level", "drift")),
    "`component` must be a single character string."
  )

  expect_error(
    autoplot(res_tempssm, component = "unknown"),
    "`component` must be one of: level, drift, season, ar1"
  )
})
