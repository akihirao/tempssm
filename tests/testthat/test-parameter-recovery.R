fit_known_exogenous_effect <- function(seed,
                                      noise_sd,
                                      add_trivial_noise = FALSE,
                                      inits = NULL) {
  set.seed(seed)

  n <- 240
  true_beta <- 2
  exo <- rep(c(-1, 0, 1, 0), length.out = n)
  temp <- 10 + true_beta * exo + stats::rnorm(n, mean = 0, sd = noise_sd)

  if (add_trivial_noise) {
    eps_noise <- .Machine$double.eps * max(1, max(abs(temp)))
    temp <- temp + eps_noise * rep(c(-1, 1), length.out = n)
  }

  temp_ts <- stats::ts(
    temp,
    start = c(2000, 1),
    frequency = 4
  )
  exo_ts <- stats::ts(
    exo,
    start = c(2000, 1),
    frequency = 4
  )
  exo_ts <- set_ts_name(exo_ts, label = "exo", quiet = TRUE)

  res <- tempssm(
    temp_data = temp_ts,
    exo_data = exo_ts,
    use_season = FALSE,
    ar_order = 1,
    inits = inits,
    maxit = 3000,
    reltol = 1e-10
  )

  list(
    true_beta = true_beta,
    res = res,
    coef = get_exo_coef(res)
  )
}


test_that("tempssm recovers a known exogenous coefficient", {
  fit <- fit_known_exogenous_effect(seed = 202405, noise_sd = 0.05)

  expect_true(fit$res$converged)
  expect_equal(fit$coef$Coefficient, fit$true_beta, tolerance = 0.02)
})


test_that("exogenous coefficient uncertainty increases with noise", {
  low_noise <- fit_known_exogenous_effect(seed = 202405, noise_sd = 0.5)
  high_noise <- fit_known_exogenous_effect(seed = 202405, noise_sd = 1.0)

  low_width <- low_noise$coef$upr - low_noise$coef$lwr
  high_width <- high_noise$coef$upr - high_noise$coef$lwr

  expect_true(low_noise$res$converged)
  expect_true(high_noise$res$converged)
  expect_gt(high_width, low_width)
})


test_that("trivial numeric noise does not materially change estimates", {
  base_fit <- fit_known_exogenous_effect(seed = 202405, noise_sd = 0.5)
  perturbed_fit <- fit_known_exogenous_effect(
    seed = 202405,
    noise_sd = 0.5,
    add_trivial_noise = TRUE
  )

  expect_true(base_fit$res$converged)
  expect_true(perturbed_fit$res$converged)
  expect_equal(
    perturbed_fit$coef[c("Coefficient", "lwr", "upr")],
    base_fit$coef[c("Coefficient", "lwr", "upr")],
    tolerance = 1e-6
  )
  expect_equal(
    as.numeric(stats::logLik(perturbed_fit$res$model)),
    as.numeric(stats::logLik(base_fit$res$model)),
    tolerance = 1e-6
  )
})


test_that("moderately different initial values give stable estimates", {
  default_fit <- fit_known_exogenous_effect(seed = 202405, noise_sd = 0.5)
  alternate_fit <- fit_known_exogenous_effect(
    seed = 202405,
    noise_sd = 0.5,
    inits = c(-12, 0.2, -0.1, -4)
  )

  expect_true(default_fit$res$converged)
  expect_true(alternate_fit$res$converged)
  expect_equal(
    alternate_fit$coef[c("Coefficient", "lwr", "upr")],
    default_fit$coef[c("Coefficient", "lwr", "upr")],
    tolerance = 1e-5
  )
  expect_equal(
    as.numeric(stats::logLik(alternate_fit$res$model)),
    as.numeric(stats::logLik(default_fit$res$model)),
    tolerance = 1e-5
  )
})
