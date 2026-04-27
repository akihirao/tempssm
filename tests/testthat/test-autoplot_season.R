test_that("autoplot_season returns a ggplot object", {

  data(niigata_sst)
  res <- tempssm(niigata_sst)

  p1 <- autoplot_season(res)
  p2 <- autoplot_season(res, ci = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})


test_that("autoplot_season checks inputs correctly", {

  data(niigata_sst)
  res <- tempssm(niigata_sst)

  expect_error(
    autoplot_season(NULL),
    "`res` must be an object returned by tempssm"
  )

  expect_error(
    autoplot_season(res, ci_level = 1.5),
    "`ci_level` must be a numeric value between 0 and 1"
  )

  expect_error(
    autoplot_season(res, ci = "yes"),
    "`ci` must be a single logical value"
  )
})

