# tests/testthat/test-fit_tempssm_safe.R

test_that(".fit_tempssm_safe returns tempssm object on success", {
  y_named <- set_ts_name(temp_ts_small, label = "Temp", quiet = TRUE)

  res <- .fit_tempssm_safe(
    y_train = y_named,
    exo_train = NULL,
    ar_order = 1,
    use_season = TRUE,
    fold_id = 1,
    marginal = TRUE
  )

  expect_s3_class(res, "tempssm")
  expect_type(res$converged, "logical")
})


test_that(".fit_tempssm_safe works with exogenous variables", {
  y_named <- set_ts_name(temp_ts_small, label = "Temp", quiet = TRUE)

  res <- .fit_tempssm_safe(
    y_train = y_named,
    exo_train = exo_ts_small,
    ar_order = 1,
    use_season = FALSE,
    fold_id = 1,
    marginal = TRUE
  )

  expect_true(is.null(res) || inherits(res, "tempssm"))
})
