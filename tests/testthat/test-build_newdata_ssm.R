# tests/testthat/test-build_newdata_ssm.R

test_that(".build_newdata_ssm works with seasonal model", {
  x <- window(exo_ts_small, end = c(2002, 12))

  exo_mat <- as.matrix(x)
  n_ahead <- NROW(exo_mat)

  pars <- rep(0.1, 6) # assuming ar=1

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    n_ahead,
    freq = 12,
    ar_order = 1,
    use_season = TRUE
  )

  expect_s3_class(model, "SSModel")
})


test_that(".build_newdata_ssm works without seasonal model", {
  x <- window(exo_ts_small, end = c(2002, 12))

  exo_mat <- as.matrix(x)
  n_ahead <- NROW(exo_mat)

  pars <- rep(0.1, 5)

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    n_ahead,
    freq = 12,
    ar_order = 1,
    use_season = FALSE
  )

  expect_s3_class(model, "SSModel")
})


test_that("number of observations matches forecast horizon", {
  x <- window(exo_ts_small, end = c(2002, 12))

  exo_mat <- as.matrix(x)
  n_ahead <- NROW(exo_mat)

  pars <- rep(0.1, 6)

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    n_ahead,
    freq = 12,
    ar_order = 1,
    use_season = TRUE
  )

  expect_identical(NROW(model$y), n_ahead)
})


test_that("observation variance is positive", {
  x <- window(exo_ts_small, end = c(2002, 12))

  exo_mat <- as.matrix(x)
  n_ahead <- NROW(exo_mat)

  pars <- rep(0.1, 6)

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    n_ahead,
    freq = 12,
    ar_order = 1,
    use_season = TRUE
  )

  expect_true(all(model$H > 0))
})


test_that("state dimension differs with seasonal option", {
  x <- window(exo_ts_small, end = c(2002, 12))

  exo_mat <- as.matrix(x)
  n_ahead <- NROW(exo_mat)

  pars_season <- rep(0.1, 6)
  pars_no <- rep(0.1, 5)

  model_season <- .build_newdata_ssm(
    pars_season,
    exo_mat,
    n_ahead,
    12,
    1,
    TRUE
  )

  model_no <- .build_newdata_ssm(
    pars_no,
    exo_mat,
    n_ahead,
    12,
    1,
    FALSE
  )

  expect_gt(ncol(model_season$T), ncol(model_no$T))
})


test_that("newdata response is univariate and unknown", {
  x <- window(exo_ts_small, end = c(2002, 12))

  exo_mat <- as.matrix(x)
  n_ahead <- NROW(exo_mat)

  pars <- rep(0.1, 6)

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    n_ahead,
    freq = 12,
    ar_order = 1,
    use_season = TRUE
  )

  expect_identical(NCOL(model$y), 1L)
  expect_true(all(is.na(model$y)))
})


test_that(".build_newdata_ssm checks exogenous horizon", {
  x <- window(exo_ts_small, end = c(2002, 12))
  exo_mat <- as.matrix(x)

  expect_error(
    .build_newdata_ssm(
      pars = rep(0.1, 6),
      exo_mat = exo_mat,
      n_ahead = NROW(exo_mat) + 1L,
      freq = 12,
      ar_order = 1,
      use_season = TRUE
    ),
    "`exo_mat` must have `n_ahead` rows."
  )
})
