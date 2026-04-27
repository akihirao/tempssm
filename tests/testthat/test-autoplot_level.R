test_that("autoplot_level returns a ggplot object", {

  data(niigata_sst)
  res <- tempssm(niigata_sst)

  p1 <- autoplot_level(res)
  p2 <- autoplot_level(res, ci = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("autoplot_level checks inputs correctly", {

  data(niigata_sst)
  res <- tempssm(niigata_sst)

  expect_error(
    autoplot_level(NULL),
    "`res` must be an object returned by tempssm"
  )

  expect_error(
    autoplot_level(res, ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1"
  )

  expect_error(
    autoplot_level(res, ci = "yes"),
    "`ci` must be a single logical value"
  )
})

