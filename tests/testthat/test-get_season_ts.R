#tests/testthat/test-get_season_ts.R

test_that("get_season_ts returns a ts object", {

  ts_obj <- get_season_ts(res_tempssm)

  expect_s3_class(ts_obj, "ts")
  expect_equal(NCOL(ts_obj), 1)
})


test_that("get_season_ts returns CI columns when ci = TRUE", {

  ts_ci <- get_season_ts(res_tempssm, ci = TRUE)

  expect_s3_class(ts_ci, "ts")
  expect_equal(colnames(ts_ci), c("season", "lwr", "upr"))
})


test_that("get_season_ts checks inputs correctly", {

  expect_error(
    get_season_ts("not a model"),
    "`res` must be a tempssm object"
  )
})
