test_that("autoplot.tempssm works for all components", {
  p1 <- autoplot(res_tempssm, component = "level")
  expect_s3_class(p1, "ggplot")

  p_all <- autoplot(res_tempssm)
  expect_true(inherits(p_all, "patchwork"))
})


test_that("autoplot works for each component", {

  components <- c("level", "drift", "season", "ar1")

  for (comp in components) {
    p <- autoplot(res_tempssm, component = comp)
    expect_s3_class(p, "ggplot")
  }
})


test_that("invalid component throws error", {

  expect_error(
    autoplot(res_tempssm, component = "invalid"),
    "must be one of"
  )
})


test_that("component must be single character", {

  expect_error(
    autoplot(res_tempssm, component = 1),
    "must be a single character"
  )

  expect_error(
    autoplot(res_tempssm, component = c("level", "drift")),
    "must be a single character"
  )
})


test_that("ci argument is passed correctly", {

  called_ci <- NULL

  mockery::stub(
    autoplot.tempssm,
    "autoplot_level",
    function(object, ci, ...) {
      called_ci <<- ci
      ggplot2::ggplot()
    }
  )

  autoplot(res_tempssm, component = "level", ci = FALSE)

  expect_false(called_ci)
})


test_that("ci_level is passed correctly", {

  called_level <- NULL

  mockery::stub(
    autoplot.tempssm,
    "autoplot_level",
    function(object, ci_level, ...) {
      called_level <<- ci_level
      ggplot2::ggplot()
    }
  )

  autoplot(res_tempssm, component = "level", ci_level = 0.9)

  expect_equal(called_level, 0.9)
})



test_that("additional arguments are passed", {

  called_args <- NULL

  mockery::stub(
    autoplot.tempssm,
    "autoplot_level",
    function(object, ..., ci) {
      called_args <<- list(...)
      ggplot2::ggplot()
    }
  )

  autoplot(res_tempssm, component = "level", some_arg = 123)

  expect_true("some_arg" %in% names(called_args))
})


