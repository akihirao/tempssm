# test-extractors.R

test_that("extract_level_ts returns ts of correct length", {

  res_temp_test <- tempssm(temp_ts_test)
  level <- extract_level_ts(res_temp_test)

  expect_s3_class(level, "ts")
  expect_equal(length(level), length(temp_ts_test))
})
