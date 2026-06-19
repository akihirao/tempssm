# tests/testthat/test-plot_temp_dev.R

test_that("plot_temp_dev returns a ggplot object", {
  temp_ts <- ts(
    rnorm(36, mean = 10),
    start = c(2000, 1),
    frequency = 12
  )

  p <- plot_temp_dev(temp_ts)

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Temperature anomalies")
  expect_identical(p$labels$x, "Time")
})


test_that("plot_temp_dev supports non-monthly seasonal frequencies", {
  temp_ts <- ts(rnorm(16), frequency = 4)

  p <- plot_temp_dev(temp_ts)

  expect_s3_class(p, "ggplot")
})


test_that("plot_temp_dev checks input class", {
  expect_error(
    plot_temp_dev(1:12),
    "`ts` must be an object of class <ts>."
  )
})
