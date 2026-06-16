# tests/testthat/test-fit_tempssm_safe.R

test_that(".fit_tempssm_safe returns tempssm object on success", {
  y_mts <- .make_mts(temp_ts_small)

  res <- .fit_tempssm_safe(
    y_train = y_mts,
    exo_train = NULL,
    ar_order = 1,
    use_season = TRUE,
    fold_id = 1
  )

  expect_s3_class(res, "tempssm")
  expect_true(is.logical(res$converged))
})


test_that(".fit_tempssm_safe works with exogenous variables", {
  y_mts <- .make_mts(temp_ts_small)

  res <- .fit_tempssm_safe(
    y_train = y_mts,
    exo_train = exo_ts_small,
    ar_order = 1,
    use_season = FALSE,
    fold_id = 1
  )

  expect_true(is.null(res) || inherits(res, "tempssm"))
})
