# test-tempssm-structure.R

test_that("tempssm object structure is consistent", {
  expect_named(
    res_tempssm,
    c(
      "model",
      "fit",
      "kfs",
      "temp_data",
      "exogenous_data",
      "ar_order",
      "use_season",
      "marginal",
      "call",
      "converged",
      "state_map"
    ),
    ignore.order = TRUE
  )

  expect_s3_class(res_tempssm$model, "SSModel")
  expect_s3_class(res_tempssm$kfs, "KFS")
  expect_s3_class(res_tempssm$temp_data, "ts")
  expect_false(res_tempssm$marginal)

  alpha <- res_tempssm$kfs$alphahat
  filtered_alpha <- res_tempssm$kfs$att
  expect_true(is.matrix(alpha))
  expect_true(is.matrix(filtered_alpha))
  expect_identical(nrow(alpha), length(res_tempssm$temp_data))
  expect_identical(nrow(filtered_alpha), length(res_tempssm$temp_data))
  expect_identical(colnames(filtered_alpha), colnames(alpha))

  state_names <- colnames(alpha)
  expect_true("level" %in% state_names)
  expect_true("slope" %in% state_names)
})
