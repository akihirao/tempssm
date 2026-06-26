#' srr standards: extended tests
#'
#' @srrstats {G5.6b} Parameter recovery across multiple random seeds is
#' implemented as an extended test. In
#' `test-extended-parameter-recovery.R`, synthetic time series with the same
#' known exogenous coefficient are generated under seeds 202401, 202402, and
#' 202403. The test checks that all fits converge and that the recovered
#' coefficients remain within a defined tolerance of the known value. Because
#' these repeated fits increase test runtime, the test is skipped by default
#' and runs only when `TEMPSSM_EXTENDED_TESTS=true` is set, as described in
#' G5.10.
#'
#' @srrstats {G5.10} Extended tests are included in the regular `testthat`
#' framework and are switched on by the `TEMPSSM_EXTENDED_TESTS=true`
#' environment variable. The file `test-extended-parameter-recovery.R` is
#' discovered by `devtools::test()` like other tests, but skips by default
#' unless that environment variable is set. When enabled, it runs additional
#' parameter-recovery checks across multiple random seeds. The same flag can
#' be set in local developer sessions or in GitHub Actions when extended
#' validation is desired without increasing the cost of every pull request.
#'
#' @srrstats {G5.12} Conditions for running extended tests are documented in
#' `CONTRIBUTING.md`. The documentation describes the
#' `TEMPSSM_EXTENDED_TESTS=true` switch, shows the command for running
#' extended tests through the standard `testthat` framework, and states that
#' the current extended tests require no downloads, large external data sets,
#' special platforms, or manual inspection of generated artefacts. It also
#' notes that the current extended test adds only a few seconds on a typical
#' development machine.
#'
test_that("extended parameter recovery is stable across seeds", {
  skip_if_not(
    identical(Sys.getenv("TEMPSSM_EXTENDED_TESTS"), "true"),
    "Set TEMPSSM_EXTENDED_TESTS=true to run extended tests."
  )

  recover_exo_coef <- function(seed) {
    set.seed(seed)

    n <- 240
    true_beta <- 2
    exo <- rep(c(-1, 0, 1, 0), length.out = n)
    temp <- 10 + true_beta * exo + stats::rnorm(n, mean = 0, sd = 0.5)

    temp_ts <- stats::ts(temp, start = c(2000, 1), frequency = 4)
    exo_ts <- stats::ts(exo, start = c(2000, 1), frequency = 4)
    exo_ts <- set_ts_name(exo_ts, label = "exo", quiet = TRUE)

    res <- tempssm(
      temp_data = temp_ts,
      exo_data = exo_ts,
      use_season = FALSE,
      ar_order = 1,
      maxit = 3000,
      reltol = 1e-10
    )

    list(
      converged = res$converged,
      estimate = get_exo_coef(res)$Coefficient
    )
  }

  fits <- lapply(c(202401, 202402, 202403), recover_exo_coef)
  estimates <- vapply(fits, `[[`, numeric(1), "estimate")
  converged <- vapply(fits, `[[`, logical(1), "converged")

  expect_true(all(converged))
  expect_equal(estimates, rep(2, length(estimates)), tolerance = 0.05)
})
