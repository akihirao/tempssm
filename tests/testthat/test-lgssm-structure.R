# test-lgssm-structure.R

test_that("ThermoSSM object structure is consistent", {

  res <- lgssm(temp_ts_test)

  expect_named(res, c("model", "fit", "kfs", "data_temp","data_exogenous","ar_order","call"), ignore.order = TRUE)

  expect_s3_class(res$model, "SSModel")
  expect_s3_class(res$kfs, "KFS")
  expect_s3_class(res$data_temp, "ts")

  alpha <- res$kfs$alphahat
  expect_true(is.matrix(alpha))
  expect_equal(nrow(alpha), length(res$data_temp))

  state_names <- colnames(alpha)
  expect_true("level" %in% state_names)
  expect_true("slope" %in% state_names)
})
