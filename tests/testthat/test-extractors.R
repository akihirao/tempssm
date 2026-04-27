# test-extractors.R

test_that("extract_level_ts returns ts of correct length", {

  level <- extract_level_ts(res_tempssm)

  expect_s3_class(level, "ts")
  expect_equal(length(level), length(temp_ts_test))
})
