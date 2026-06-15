# tests/testthat/test-get_ar1_ts.R

test_that("get_ar1_ts basic structure", {
  ts_obj <- get_ar1_ts(res_tempssm)

  check_ts_basic(ts_obj, res_tempssm)
  check_ts_univariate(ts_obj)
})

test_that("get_ar1_ts CI structure", {
  ts_ci <- get_ar1_ts(res_tempssm, ci = TRUE)

  check_ts_ci_structure(ts_ci, "ar1")
})

test_that("get_ar1_ts CI values", {
  ts_ci <- get_ar1_ts(res_tempssm, ci = TRUE)

  ci_obj <- stats::confint(res_tempssm$kfs)

  check_ts_ci_values(ts_ci, ci_obj, "arima1")
})
