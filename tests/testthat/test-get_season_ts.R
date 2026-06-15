# tests/testthat/test-get_season_ts.R

# ---- basic & CI --------------------------------------------------------

test_that("get_season_ts basic and CI behaviour", {

  # basic
  ts_obj <- get_season_ts(res_tempssm)
  check_ts_basic(ts_obj, res_tempssm)
  check_ts_univariate(ts_obj)

  # CI structure
  ts_ci <- get_season_ts(res_tempssm, ci = TRUE)
  check_ts_ci_structure(ts_ci, "season")

  # CI values
  ci_obj <- stats::confint(res_tempssm$kfs)
  check_ts_ci_values(ts_ci, ci_obj, "sea_dummy1")
})


# ---- seasonal-specific logic -------------------------------------------

test_that("get_season_ts requires seasonal model", {

  bad_res <- res_tempssm
  bad_res$use_season <- FALSE

  expect_error(
    get_season_ts(bad_res),
    "Seasonal component is not included"
  )
})


# ---- smoothing structure ------------------------------------------------

test_that("get_season_ts validates smoothing components", {

  # alphahat NULL
  bad_res1 <- res_tempssm
  bad_res1$kfs$alphahat <- NULL

  expect_error(
    get_season_ts(bad_res1),
    "Seasonal component not found"
  )

  # missing column
  bad_res2 <- res_tempssm
  bad_res2$kfs$alphahat <- matrix(rnorm(10), ncol = 1)
  colnames(bad_res2$kfs$alphahat) <- "level"

  expect_error(
    get_season_ts(bad_res2),
    "Seasonal component not found"
  )
})


# ---- CI structure errors -----------------------------------------------

test_that("get_season_ts validates CI content", {

  bad_res <- res_tempssm

  testthat::local_mocked_bindings(
    confint = function(...) list(level = matrix(0)),
    .package = "stats"
  )

  expect_error(
    get_season_ts(bad_res, ci = TRUE),
    "Seasonal component not found in confidence intervals"
  )
})


# ---- input validation ---------------------------------------------------

test_that("get_season_ts input validation", {
  check_invalid_input(get_season_ts, res_tempssm)
})
