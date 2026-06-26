#' srr standards: parameter recovery and noise susceptibility
#'
#' @srrstats {G5.6} Parameter recovery tests are included for exogenous
#' regression effects, which are directly interpretable user-facing
#' parameters in `tempssm` models. In `test-parameter-recovery.R`, a synthetic
#' regular `ts` series is generated from a known exogenous coefficient of 2
#' with small Gaussian noise. The test fits `tempssm()` with the corresponding
#' exogenous `ts` input and checks that `get_exo_coef()` recovers the known
#' coefficient within tolerance. Variance, Kalman filtering, and smoothing
#' internals are delegated to `KFAS` and are checked separately through direct
#' backend comparisons.
#'
#' @srrstats {G5.6a} Parameter recovery tests use explicit numerical
#' tolerances rather than exact equality. In `test-parameter-recovery.R`, the
#' known exogenous coefficient is 2 and the recovered coefficient from
#' `get_exo_coef()` is checked with `expect_equal(..., tolerance = 0.02)`.
#' This reflects the fact that `tempssm()` estimates state-space parameters by
#' numerical optimization and Kalman filtering, so recovered values are
#' expected to agree within a defined tolerance rather than match exactly.
#'
#' @srrstats {G5.7} Algorithm performance is tested against changing data
#' properties in `test-parameter-recovery.R`. Using the same known exogenous
#' coefficient and fixed random seed, the test fits one synthetic series with
#' moderate Gaussian noise and another with higher Gaussian noise. It checks
#' that both models converge and that the confidence interval width for the
#' recovered exogenous coefficient is larger in the higher-noise case, as
#' expected for a less informative time series.
#'
#' @srrstats {G5.9} Noise susceptibility is tested with synthetic time-series
#' data in `test-parameter-recovery.R`. The tests generate a known exogenous
#' coefficient under fixed random seeds, verify coefficient recovery under
#' small Gaussian noise, and compare moderate- and higher-noise versions of
#' the same design. The higher-noise series is expected to yield a wider
#' confidence interval for the recovered exogenous coefficient while both
#' fits converge, demonstrating expected stochastic behaviour as observation
#' noise changes.
#'
#' @srrstats {G5.9a} The test suite checks that adding deterministic
#' `.Machine$double.eps`-scale noise to a synthetic temperature series does
#' not materially change fitted results. In `test-parameter-recovery.R`, the
#' same known exogenous-effect series is fitted with and without a tiny
#' alternating perturbation. The test verifies that both models converge and
#' that the recovered exogenous coefficient, confidence limits, and fitted
#' log-likelihood remain equal within numerical tolerance.
#'
#' @srrstats {G5.9b} Stability under different initial conditions is tested
#' in `test-parameter-recovery.R`. The same synthetic time series with a
#' known exogenous coefficient is fitted once with default initial values and
#' once with a moderately different user-supplied `inits` vector. The test
#' verifies that both fits converge and that the recovered exogenous
#' coefficient, confidence limits, and fitted log-likelihood remain equal
#' within numerical tolerance. The fitting algorithm itself is deterministic
#' for fixed data and initial values, so random seeds affect only synthetic
#' data generation in these tests.
#'
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
