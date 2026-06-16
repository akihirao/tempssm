# tests/testthat/test-build_newdata_ssm.R

test_that(".build_newdata_ssm works with seasonal model", {
  y <- window(temp_ts_small, end = c(2002, 12))
  x <- window(exo_ts_small, end = c(2002, 12))

  y_mts <- .make_mts(y)

  exo_mat <- as.matrix(x)
  data <- cbind(y_mts, x)

  pars <- rep(0.1, 6) # assuming ar=1

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    data,
    freq = 12,
    ar_order = 1,
    use_season = TRUE
  )

  expect_s3_class(model, "SSModel")
})


test_that(".build_newdata_ssm works without seasonal model", {
  y <- window(temp_ts_small, end = c(2002, 12))
  x <- window(exo_ts_small, end = c(2002, 12))

  y_mts <- .make_mts(y)
  exo_mat <- as.matrix(x)
  data <- cbind(y_mts, x)

  pars <- rep(0.1, 5)

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    data,
    freq = 12,
    ar_order = 1,
    use_season = FALSE
  )

  expect_s3_class(model, "SSModel")
})


test_that("number of observations matches input data", {
  y <- window(temp_ts_small, end = c(2002, 12))
  x <- window(exo_ts_small, end = c(2002, 12))

  y_mts <- .make_mts(y)
  exo_mat <- as.matrix(x)
  data <- cbind(y_mts, x)

  pars <- rep(0.1, 6)

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    data,
    freq = 12,
    ar_order = 1,
    use_season = TRUE
  )

  expect_equal(NROW(model$y), NROW(data))
})


test_that("observation variance is positive", {
  y <- window(temp_ts_small, end = c(2002, 12))
  x <- window(exo_ts_small, end = c(2002, 12))

  y_mts <- .make_mts(y)
  exo_mat <- as.matrix(x)
  data <- cbind(y_mts, x)

  pars <- rep(0.1, 6)

  model <- .build_newdata_ssm(
    pars,
    exo_mat,
    data,
    freq = 12,
    ar_order = 1,
    use_season = TRUE
  )

  expect_true(all(model$H > 0))
})


test_that("state dimension differs with seasonal option", {
  y <- window(temp_ts_small, end = c(2002, 12))
  x <- window(exo_ts_small, end = c(2002, 12))

  y_mts <- .make_mts(y)
  exo_mat <- as.matrix(x)
  data <- cbind(y_mts, x)

  pars_season <- rep(0.1, 6)
  pars_no <- rep(0.1, 5)

  model_season <- .build_newdata_ssm(
    pars_season,
    exo_mat,
    data,
    12,
    1,
    TRUE
  )

  model_no <- .build_newdata_ssm(
    pars_no,
    exo_mat,
    data,
    12,
    1,
    FALSE
  )

  expect_true(ncol(model_season$T) > ncol(model_no$T))
})


test_that("handles NA data input", {
  y <- window(temp_ts_small, end = c(2002, 12))
  y[] <- NA

  x <- window(exo_ts_small, end = c(2002, 12))

  y_mts <- .make_mts(y)
  exo_mat <- as.matrix(x)
  data <- cbind(y_mts, x)

  pars <- rep(0.1, 6)

  expect_no_error(
    .build_newdata_ssm(
      pars,
      exo_mat,
      data,
      freq = 12,
      ar_order = 1,
      use_season = TRUE
    )
  )
})
