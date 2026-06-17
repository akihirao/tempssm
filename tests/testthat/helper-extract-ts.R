# tests/testthat/helper-extract-ts.R

# ---- basic structure checks --------------------------------------------

check_ts_basic <- function(ts_obj, res) {
  testthat::expect_s3_class(ts_obj, "ts")
  testthat::expect_identical(start(ts_obj), start(res$temp_data))
  testthat::expect_identical(frequency(ts_obj), frequency(res$temp_data))
}

check_ts_univariate <- function(ts_obj) {
  expect_equal(NCOL(ts_obj), 1)
}


# ---- CI checks ---------------------------------------------------------


check_ts_ci_structure <- function(ts_ci, value_name) {
  testthat::expect_s3_class(ts_ci, "ts")
  testthat::expect_equal(NCOL(ts_ci), 3)
  testthat::expect_named(
    as.data.frame(ts_ci),
    c(value_name, "lwr", "upr")
  )
}


check_ts_ci_values <- function(ts_ci, ci_obj, key, scale = 1) {
  testthat::expect_equal(ts_ci[, "lwr"], ci_obj[[key]][, "lwr"] * scale)
  testthat::expect_equal(ts_ci[, "upr"], ci_obj[[key]][, "upr"] * scale)
}

# ---- error helper ------------------------------------------------------

check_invalid_input <- function(fun, res) {
  testthat::expect_error(
    fun("not a model"),
    "`res` must be an object of class"
  )

  testthat::expect_error(
    fun(res, ci = TRUE, ci_level = 1.2),
    "`ci_level` must be a numeric value between 0 and 1"
  )
}

# ---- smoothing component missing ---------------------------------------

check_missing_component <- function(fun, res, colname, msg) {
  bad_res <- res
  bad_res$kfs$alphahat <- matrix(rnorm(10), ncol = 1)
  colnames(bad_res$kfs$alphahat) <- "dummy"

  testthat::expect_error(fun(bad_res), msg)
}
