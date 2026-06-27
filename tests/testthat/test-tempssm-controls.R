test_that(".prepare_tempssm_numeric_control handles defaults and validation", {
  expect_identical(
    .prepare_tempssm_numeric_control(NULL, "maxit", 5000),
    5000
  )
  expect_identical(
    .prepare_tempssm_numeric_control(1000, "maxit", 5000),
    1000
  )
  expect_error(
    .prepare_tempssm_numeric_control(c(100, 200), "maxit", 5000),
    "maxit.*length one"
  )
  expect_error(
    .prepare_tempssm_numeric_control(Inf, "maxit", 5000),
    "maxit.*finite numeric scalar"
  )
})


test_that(".prepare_tempssm_controls supplies defaults", {
  controls <- .prepare_tempssm_controls(
    ar_order = 1,
    use_season = TRUE,
    maxit = NULL,
    reltol = NULL,
    na_action = "inform"
  )

  expect_identical(
    controls,
    list(
      ar_order = 1,
      use_season = TRUE,
      maxit = 5000,
      reltol = 1e-16,
      na_action = "inform"
    )
  )
})


test_that(".prepare_tempssm_controls preserves valid supplied controls", {
  controls <- .prepare_tempssm_controls(
    ar_order = 3,
    use_season = FALSE,
    maxit = 1000,
    reltol = 1e-10,
    na_action = "allow"
  )

  expect_identical(controls$ar_order, 3)
  expect_false(controls$use_season)
  expect_identical(controls$maxit, 1000)
  expect_identical(controls$reltol, 1e-10)
  expect_identical(controls$na_action, "allow")
})


test_that(".prepare_tempssm_controls rejects invalid controls", {
  expect_error(
    .prepare_tempssm_controls(0, TRUE, NULL, NULL, "inform"),
    "ar_order.*integer >= 1"
  )
  expect_error(
    .prepare_tempssm_controls(1, NA, NULL, NULL, "inform"),
    "use_season.*logical scalar"
  )
  expect_error(
    .prepare_tempssm_controls(1, TRUE, Inf, NULL, "inform"),
    "maxit.*finite numeric scalar"
  )
  expect_error(
    .prepare_tempssm_controls(1, TRUE, NULL, NA_real_, "inform"),
    "reltol.*finite numeric scalar"
  )
  expect_error(
    .prepare_tempssm_controls(1, TRUE, NULL, NULL, "drop"),
    "should be one of"
  )
})


test_that(".prepare_tempssm_controls warns for high AR orders", {
  expect_warning(
    controls <- .prepare_tempssm_controls(
      ar_order = 5,
      use_season = TRUE,
      maxit = NULL,
      reltol = NULL,
      na_action = "inform"
    ),
    "ar_order.*greater than 4"
  )

  expect_identical(controls$ar_order, 5)
})
