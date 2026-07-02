make_future_exogenous <- function(values,
                                   start = c(2010, 1),
                                   frequency = 12) {
  ts(
    matrix(
      values,
      ncol = 1,
      dimnames = list(
        NULL,
        colnames(res_tempssm_exo$exogenous_data)
      )
    ),
    start = start,
    frequency = frequency
  )
}


test_that("predict.tempssm forecasts the next time point by default", {
  pred <- predict(res_tempssm)

  expect_s3_class(pred, "ts")
  expect_length(pred, 1L)
  expect_identical(
    as.numeric(stats::time(pred)[1L]),
    as.numeric(tail(stats::time(res_tempssm$temp_data), 1L)) +
      1 / stats::frequency(res_tempssm$temp_data)
  )
})


test_that("predict.tempssm agrees with the KFAS model forecast", {
  pred <- predict(res_tempssm, n.ahead = 4L)
  expected <- stats::predict(res_tempssm$model, n.ahead = 4L)

  expect_identical(pred, expected)
})


test_that("predict.tempssm returns requested prediction intervals", {
  pred <- predict(
    res_tempssm,
    n.ahead = 3L,
    interval = "prediction",
    level = 0.9
  )

  expect_s3_class(pred, "mts")
  expect_identical(colnames(pred), c("fit", "lwr", "upr"))
  expect_identical(NROW(pred), 3L)
  expect_true(all(pred[, "lwr"] <= pred[, "fit"]))
  expect_true(all(pred[, "fit"] <= pred[, "upr"]))
})


test_that("predict.tempssm forecasts with future exogenous values", {
  future_exogenous <- make_future_exogenous(c(0.1, 0.2, 0.3))

  pred <- predict(
    res_tempssm_exo,
    new_exo_data = future_exogenous
  )
  expected <- .predict_with_exo(res_tempssm_exo, future_exogenous)

  expect_s3_class(pred, "ts")
  expect_length(pred, 3L)
  expect_identical(pred, expected)
  expect_identical(stats::start(pred), c(2010, 1))
})


test_that("predict.tempssm returns intervals for exogenous forecasts", {
  future_exogenous <- make_future_exogenous(c(0.1, 0.2))

  pred <- predict(
    res_tempssm_exo,
    new_exo_data = future_exogenous,
    interval = "prediction"
  )

  expect_s3_class(pred, "mts")
  expect_identical(colnames(pred), c("fit", "lwr", "upr"))
  expect_identical(NROW(pred), 2L)
})


test_that("predict.tempssm carries exogenous values forward explicitly", {
  last_values <- as.numeric(
    res_tempssm_exo$exogenous_data[
      NROW(res_tempssm_exo$exogenous_data),
      ,
      drop = TRUE
    ]
  )
  explicit_future <- make_future_exogenous(last_values)
  expected <- predict(
    res_tempssm_exo,
    new_exo_data = explicit_future
  )

  withr::local_options(list(tempssm.verbosity = "inform"))
  expect_message(
    pred <- predict(res_tempssm_exo, exo_strategy = "last"),
    "persistence assumption"
  )

  expect_identical(pred, expected)
  expect_length(pred, 1L)
})


test_that("last exogenous strategy is restricted to its stated scope", {
  future_exogenous <- make_future_exogenous(0.1)

  expect_error(
    predict(
      res_tempssm_exo,
      new_exo_data = future_exogenous,
      exo_strategy = "last"
    ),
    "new_exo_data.*NULL"
  )
  expect_error(
    predict(
      res_tempssm_exo,
      n.ahead = 2L,
      exo_strategy = "last"
    ),
    "limited to one-step"
  )
  expect_error(
    predict(res_tempssm, exo_strategy = "last"),
    "only be used with an exogenous model"
  )
})


test_that("last exogenous strategy preserves multiple covariates", {
  fitted_exogenous <- ts(
    matrix(
      1:8,
      ncol = 2,
      dimnames = list(NULL, c("x1", "x2"))
    ),
    start = c(2000, 1),
    frequency = 4
  )
  object <- structure(
    list(exogenous_data = fitted_exogenous),
    class = "tempssm"
  )

  future_exogenous <- .make_last_exogenous_forecast(object)

  expect_s3_class(future_exogenous, "mts")
  expect_identical(as.numeric(future_exogenous), c(4, 8))
  expect_identical(colnames(future_exogenous), c("x1", "x2"))
  expect_identical(stats::start(future_exogenous), c(2001, 1))
})


test_that("predict.tempssm requires covariates only for exogenous models", {
  expect_error(
    predict(res_tempssm_exo),
    "new_exo_data.*required"
  )
  expect_error(
    predict(
      res_tempssm,
      new_exo_data = ts(1, start = c(2010, 1), frequency = 12)
    ),
    "new_exo_data.*NULL"
  )
})


test_that("predict.tempssm validates future exogenous structure", {
  valid <- make_future_exogenous(c(0.1, 0.2))

  expect_error(
    predict(res_tempssm_exo, n.ahead = 1L, new_exo_data = valid),
    "number of rows.*n.ahead"
  )
  expect_error(
    predict(res_tempssm_exo, new_exo_data = unclass(valid)),
    "new_exo_data.*ts"
  )

  missing_value <- valid
  missing_value[1L, 1L] <- NA_real_
  expect_error(
    predict(res_tempssm_exo, new_exo_data = missing_value),
    "missing or non-finite"
  )

  wrong_name <- valid
  colnames(wrong_name) <- "wrong"
  expect_error(
    predict(res_tempssm_exo, new_exo_data = wrong_name),
    "Column names and order"
  )

  wrong_frequency <- make_future_exogenous(
    c(0.1, 0.2),
    frequency = 4
  )
  expect_error(
    predict(res_tempssm_exo, new_exo_data = wrong_frequency),
    "Frequency"
  )

  wrong_start <- stats::window(valid, start = c(2010, 2))
  expect_error(
    predict(res_tempssm_exo, new_exo_data = wrong_start),
    "first time point after"
  )
})


test_that("predict.tempssm validates the fitted object", {
  expect_error(
    predict(structure(list(), class = "tempssm")),
    "requires a converged"
  )
})


test_that("predict.tempssm validates forecast controls", {
  invalid_horizons <- list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "1")
  for (horizon in invalid_horizons) {
    expect_error(
      predict(res_tempssm, n.ahead = horizon),
      "n.ahead.*single positive integer"
    )
  }

  invalid_levels <- list(0, 1, NA_real_, Inf, c(0.8, 0.9), "0.95")
  for (confidence_level in invalid_levels) {
    expect_error(
      predict(res_tempssm, level = confidence_level),
      "level.*between 0 and 1"
    )
  }

  expect_error(
    predict(res_tempssm, interval = "unsupported"),
    "arg.*should be one of"
  )
})
