test_that("predict.tempssm forecasts the next time point by default", {
  pred <- predict(res_tempssm)

  expect_s3_class(pred, "ts")
  expect_length(pred, 1L)
  expect_equal(
    as.numeric(stats::time(pred)[1L]),
    as.numeric(tail(stats::time(res_tempssm$temp_data), 1L)) +
      1 / stats::frequency(res_tempssm$temp_data)
  )
})


test_that("predict.tempssm agrees with the KFAS model forecast", {
  pred <- predict(res_tempssm, n.ahead = 4L)
  expected <- stats::predict(res_tempssm$model, n.ahead = 4L)

  expect_equal(pred, expected)
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


test_that("predict.tempssm rejects unsupported exogenous forecasts", {
  expect_error(
    predict(res_tempssm_exo),
    "exogenous variables is not yet supported"
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
