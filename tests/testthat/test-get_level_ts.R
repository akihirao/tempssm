# tests/testthat/test-get_level_ts.R

test_that("get_level_ts returns a ts object", {
  ts_obj <- get_level_ts(res_tempssm)

  expect_s3_class(ts_obj, "ts")
  expect_identical(NCOL(ts_obj), 1L)
})


test_that("get_level_ts returns CI columns when ci = TRUE", {
  ts_ci <- get_level_ts(res_tempssm, ci = TRUE)

  expect_s3_class(ts_ci, "ts")
  expect_identical(colnames(ts_ci), c("level", "lwr", "upr"))
})


test_that("get_level_ts checks inputs correctly", {
  expect_error(
    get_level_ts(NULL),
    "`res` must be an object of class 'tempssm'."
  )

  expect_error(
    get_level_ts(res_tempssm, ci = TRUE, ci_level = 1.2),
    "`ci_level` must be a numeric value between 0 and 1."
  )

  expect_error(
    get_level_ts(res_tempssm, estimate = "predicted"),
    "arg.*should be one of"
  )
})


test_that("filtered estimates mask the diffuse phase", {
  filtered <- get_level_ts(res_tempssm, estimate = "filtered")

  expect_s3_class(filtered, "ts")
  expect_identical(NCOL(filtered), 1L)
  expect_true(all(is.na(filtered[seq_len(res_tempssm$kfs$d)])))
})


test_that("filtered confidence intervals mask the diffuse phase", {
  filtered <- get_level_ts(
    res_tempssm,
    ci = TRUE,
    estimate = "filtered"
  )
  diffuse_idx <- seq_len(res_tempssm$kfs$d)

  expect_s3_class(filtered, "ts")
  expect_identical(colnames(filtered), c("level", "lwr", "upr"))
  expect_true(all(is.na(filtered[diffuse_idx, ])))
})


test_that("filtered estimates validate the diffuse phase", {
  bad_res <- res_tempssm
  bad_res$kfs$d <- NULL

  expect_error(
    get_level_ts(bad_res, estimate = "filtered"),
    "diffuse filtering phase is unavailable or invalid"
  )

  no_diffuse <- res_tempssm
  no_diffuse$kfs$d <- 0L
  filtered <- get_level_ts(no_diffuse, estimate = "filtered")
  expect_identical(
    as.numeric(filtered),
    as.numeric(no_diffuse$kfs$att[, "level"])
  )

  filtered_ci <- get_level_ts(
    no_diffuse,
    ci = TRUE,
    estimate = "filtered"
  )
  expect_false(anyNA(filtered_ci))
})


test_that("get_level_ts validates filtered state availability", {
  bad_res <- res_tempssm
  bad_res$kfs$att <- NULL

  expect_error(
    get_level_ts(bad_res, estimate = "filtered"),
    "Level component not found in the filtering results"
  )
})


test_that("get_level_ts validates filtered covariance results", {
  missing_covariance <- res_tempssm
  missing_covariance$kfs$Ptt <- NULL
  expect_error(
    get_level_ts(
      missing_covariance,
      ci = TRUE,
      estimate = "filtered"
    ),
    "Filtered state covariance results are unavailable or invalid"
  )

  negative_variance <- res_tempssm
  state_idx <- match("level", colnames(negative_variance$kfs$att))
  time_idx <- negative_variance$kfs$d + 1L
  negative_variance$kfs$Ptt[state_idx, state_idx, time_idx] <- -1
  expect_error(
    get_level_ts(
      negative_variance,
      ci = TRUE,
      estimate = "filtered"
    ),
    "Filtered state variances must not be negative"
  )

  nonfinite_variance <- res_tempssm
  state_idx <- match("level", colnames(nonfinite_variance$kfs$att))
  time_idx <- nonfinite_variance$kfs$d + 1L
  nonfinite_variance$kfs$Ptt[state_idx, state_idx, time_idx] <- Inf
  expect_error(
    get_level_ts(
      nonfinite_variance,
      ci = TRUE,
      estimate = "filtered"
    ),
    "Filtered state variances must be finite"
  )

  rounding_variance <- res_tempssm
  state_idx <- match("level", colnames(rounding_variance$kfs$att))
  time_idx <- rounding_variance$kfs$d + 1L
  rounding_variance$kfs$Ptt[state_idx, state_idx, time_idx] <-
    -.Machine$double.eps
  rounded_ci <- get_level_ts(
    rounding_variance,
    ci = TRUE,
    estimate = "filtered"
  )
  expect_identical(
    unname(rounded_ci[time_idx, "lwr"]),
    unname(rounded_ci[time_idx, "level"])
  )
  expect_identical(
    unname(rounded_ci[time_idx, "upr"]),
    unname(rounded_ci[time_idx, "level"])
  )
})


test_that("get_level_ts preserves time attributes", {
  ts_obj <- get_level_ts(res_tempssm)

  expect_identical(start(ts_obj), start(res_tempssm$temp_data))
  expect_identical(frequency(ts_obj), frequency(res_tempssm$temp_data))
})


test_that("get_level_ts returns correct number of columns with CI", {
  ts_ci <- get_level_ts(res_tempssm, ci = TRUE)

  expect_identical(NCOL(ts_ci), 3L)
})


test_that("CI output preserves time attributes", {
  ts_ci <- get_level_ts(res_tempssm, ci = TRUE)

  expect_identical(start(ts_ci), start(res_tempssm$temp_data))
  expect_identical(frequency(ts_ci), frequency(res_tempssm$temp_data))
})


test_that("errors when level component is missing", {
  bad_res <- res_tempssm
  bad_res$kfs$alphahat <- NULL

  expect_error(
    get_level_ts(bad_res),
    "Level component not found"
  )
})


test_that("errors when level column is missing in alphahat", {
  bad_res <- res_tempssm
  bad_res$kfs$alphahat <- matrix(rnorm(10), ncol = 1)
  colnames(bad_res$kfs$alphahat) <- "trend" # no existence of "level"
  expect_error(
    get_level_ts(bad_res),
    "Level component not found"
  )
})


test_that("errors when level missing in confidence intervals", {
  bad_res <- res_tempssm

  original_confint <- stats::confint

  testthat::local_mocked_bindings(
    confint = function(...) list(other = matrix(0)),
    .package = "stats"
  )

  expect_error(
    get_level_ts(bad_res, ci = TRUE),
    "Level component not found in confidence intervals"
  )
})


test_that("ci = FALSE returns univariate ts", {
  ts_obj <- get_level_ts(res_tempssm, ci = FALSE)

  expect_identical(NCOL(ts_obj), 1L)
})


test_that("ci = FALSE returns univariate ts", {
  ts_obj <- get_level_ts(res_tempssm, ci = FALSE)

  expect_identical(NCOL(ts_obj), 1L)
})


test_that("ci_level boundary values", {
  expect_error(
    get_level_ts(res_tempssm, ci = TRUE, ci_level = 0),
    "ci_level.*between 0 and 1"
  )
  expect_error(
    get_level_ts(res_tempssm, ci = TRUE, ci_level = 1),
    "ci_level.*between 0 and 1"
  )
})
