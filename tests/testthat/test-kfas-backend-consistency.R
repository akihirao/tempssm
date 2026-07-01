test_that("tempssm logLik method matches the underlying KFAS model", {
  tempssm_loglik <- logLik(res_tempssm)
  kfas_loglik <- stats::logLik(res_tempssm$model)

  expect_identical(as.numeric(tempssm_loglik), as.numeric(kfas_loglik))
  expect_identical(attr(tempssm_loglik, "nobs"), length(res_tempssm$temp_data))
  expect_identical(
    attr(tempssm_loglik, "df"),
    length(res_tempssm$fit$optim.out$par)
  )
})


test_that("marginal logLik matches the underlying KFAS model", {
  tempssm_loglik <- logLik(res_tempssm, marginal = TRUE)
  kfas_loglik <- stats::logLik(res_tempssm$model, marginal = TRUE)

  expect_identical(as.numeric(tempssm_loglik), as.numeric(kfas_loglik))
  expect_true(attr(tempssm_loglik, "marginal"))
})


test_that("tempssm AIC uses the KFAS log-likelihood with tempssm df", {
  kfas_loglik <- as.numeric(stats::logLik(res_tempssm$model))
  expected_aic <- -2 * kfas_loglik + 2 * length(res_tempssm$fit$optim.out$par)

  expect_identical(AIC(res_tempssm), expected_aic)
  expect_identical(get_aic(res_tempssm), expected_aic)
})


test_that("component accessors match KFAS smoothed states", {
  level_ts <- get_level_ts(res_tempssm, ci = FALSE)
  drift_ts <- get_drift_ts(res_tempssm, ci = FALSE)
  season_ts <- get_season_ts(res_tempssm, ci = FALSE)
  ar1_ts <- get_ar1_ts(res_tempssm, ci = FALSE)
  kfas_states <- res_tempssm$kfs$alphahat
  freq <- stats::frequency(res_tempssm$temp_data)

  expect_identical(as.numeric(level_ts), as.numeric(kfas_states[, "level"]))
  expect_identical(
    as.numeric(drift_ts),
    as.numeric(kfas_states[, "slope"] * freq)
  )
  expect_identical(
    as.numeric(season_ts),
    as.numeric(kfas_states[, "sea_dummy1"])
  )
  expect_identical(as.numeric(ar1_ts), as.numeric(kfas_states[, "arima1"]))
})


test_that("component confidence intervals match KFAS confint output", {
  ci_obj <- stats::confint(res_tempssm$kfs, level = 0.95)
  freq <- stats::frequency(res_tempssm$temp_data)

  level_ci <- get_level_ts(res_tempssm, ci = TRUE, ci_level = 0.95)
  drift_ci <- get_drift_ts(res_tempssm, ci = TRUE, ci_level = 0.95)
  season_ci <- get_season_ts(res_tempssm, ci = TRUE, ci_level = 0.95)
  ar1_ci <- get_ar1_ts(res_tempssm, ci = TRUE, ci_level = 0.95)

  expect_identical(
    as.numeric(level_ci[, "lwr"]),
    as.numeric(ci_obj$level[, "lwr"])
  )
  expect_identical(
    as.numeric(level_ci[, "upr"]),
    as.numeric(ci_obj$level[, "upr"])
  )
  expect_identical(
    as.numeric(drift_ci[, "lwr"]),
    as.numeric(ci_obj$slope[, "lwr"] * freq)
  )
  expect_identical(
    as.numeric(drift_ci[, "upr"]),
    as.numeric(ci_obj$slope[, "upr"] * freq)
  )
  expect_identical(
    as.numeric(season_ci[, "lwr"]),
    as.numeric(ci_obj$sea_dummy1[, "lwr"])
  )
  expect_identical(
    as.numeric(season_ci[, "upr"]),
    as.numeric(ci_obj$sea_dummy1[, "upr"])
  )
  expect_identical(
    as.numeric(ar1_ci[, "lwr"]),
    as.numeric(ci_obj$arima1[, "lwr"])
  )
  expect_identical(
    as.numeric(ar1_ci[, "upr"]),
    as.numeric(ci_obj$arima1[, "upr"])
  )
})
