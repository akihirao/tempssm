# test-tempssm-structure.R

test_that("tempssm object structure is consistent", {

  expect_named(
    res_tempssm,
    c("model",
      "fit",
      "kfs",
      "data_temp",
      "data_exogenous",
      "ar_order",
      "use_season",
      "call",
      "converged",
      "state_map"
      ), ignore.order = TRUE)

  expect_s3_class(res_tempssm$model, "SSModel")
  expect_s3_class(res_tempssm$kfs, "KFS")
  expect_s3_class(res_tempssm$data_temp, "ts")

  alpha <- res_tempssm$kfs$alphahat
  expect_true(is.matrix(alpha))
  expect_equal(nrow(alpha), length(res_tempssm$data_temp))

  state_names <- colnames(alpha)
  expect_true("level" %in% state_names)
  expect_true("slope" %in% state_names)
})
