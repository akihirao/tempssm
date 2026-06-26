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
