make_exo_accessor_fixture <- function() {
  alpha_hat <- cbind(
    x1 = rep(0.25, 3),
    x2 = rep(-0.5, 3),
    level = c(10, 11, 12)
  )
  res <- list(
    converged = TRUE,
    state_map = list(exogenous = c("x1", "x2")),
    kfs = list(alphahat = alpha_hat)
  )
  class(res) <- "tempssm"

  ci <- list(
    exo_matx1 = cbind(
      lwr = rep(0.1, 3),
      upr = rep(0.4, 3)
    ),
    exo_matx2 = cbind(
      lwr = rep(-0.8, 3),
      upr = rep(-0.2, 3)
    ),
    level = cbind(
      lwr = c(9, 10, 11),
      upr = c(11, 12, 13)
    )
  )

  list(res = res, ci = ci)
}


test_that("get_exo_coef maps multiple exogenous states by KFAS position", {
  fixture <- make_exo_accessor_fixture()
  observed <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    confint = function(object, level, ...) {
      observed$object <- object
      observed$level <- level
      fixture$ci
    },
    .package = "stats"
  )

  result <- get_exo_coef(fixture$res, ci_level = 0.9)

  expect_identical(
    result,
    data.frame(
      Variable = c("x1", "x2"),
      Coefficient = c(0.25, -0.5),
      lwr = c(0.1, -0.8),
      upr = c(0.4, -0.2)
    )
  )
  expect_identical(observed$object, fixture$res$kfs)
  expect_identical(observed$level, 0.9)
})


test_that("get_exo_coef rejects unresolved exogenous states", {
  fixture <- make_exo_accessor_fixture()
  fixture$res$state_map$exogenous <- c("x1", "missing")

  expect_error(
    get_exo_coef(fixture$res),
    "states listed in `state_map` were not found"
  )
})


test_that("get_exo_coef rejects invalid confidence levels", {
  fixture <- make_exo_accessor_fixture()

  expect_error(
    get_exo_coef(fixture$res, ci_level = c(0.8, 0.9)),
    "ci_level"
  )
  expect_error(
    get_exo_coef(fixture$res, ci_level = NA_real_),
    "ci_level"
  )
  expect_error(
    get_exo_coef(fixture$res, ci_level = NaN),
    "ci_level"
  )
  expect_error(
    get_exo_coef(fixture$res, ci_level = Inf),
    "ci_level"
  )
})


test_that("get_exo_coef rejects incomplete smoothing structures", {
  fixture <- make_exo_accessor_fixture()
  missing_kfs <- fixture$res
  missing_kfs$kfs <- NULL
  missing_states <- fixture$res
  missing_states$kfs$alphahat <- NULL

  expect_error(
    get_exo_coef(missing_kfs),
    "Smoothing results for exogenous states are not available"
  )
  expect_error(
    get_exo_coef(missing_states),
    "Smoothing results for exogenous states are not available"
  )
})


test_that("get_exo_coef rejects malformed confidence intervals", {
  fixture <- make_exo_accessor_fixture()
  malformed_ci <- fixture$ci
  colnames(malformed_ci[[2]]) <- c("lower", "upper")
  testthat::local_mocked_bindings(
    confint = function(...) malformed_ci,
    .package = "stats"
  )

  expect_error(
    get_exo_coef(fixture$res),
    "Confidence intervals for all exogenous states are not available"
  )
})


test_that("get_exo_coef rejects missing confidence interval elements", {
  fixture <- make_exo_accessor_fixture()
  testthat::local_mocked_bindings(
    confint = function(...) fixture$ci[1],
    .package = "stats"
  )

  expect_error(
    get_exo_coef(fixture$res),
    "Confidence intervals for all exogenous states are not available"
  )
})


test_that("get_exo_coef rejects non-numeric confidence intervals", {
  fixture <- make_exo_accessor_fixture()
  invalid_ci <- fixture$ci
  invalid_ci[[1]][] <- as.character(invalid_ci[[1]])
  testthat::local_mocked_bindings(
    confint = function(...) invalid_ci,
    .package = "stats"
  )

  expect_error(
    get_exo_coef(fixture$res),
    "Confidence intervals for all exogenous states must be numeric"
  )
})


test_that("common accessor checks reject invalid ci flags", {
  expect_error(
    get_level_ts(res_tempssm, ci = NA),
    "ci.*non-missing logical scalar"
  )
  expect_error(
    get_level_ts(res_tempssm, ci = c(TRUE, FALSE)),
    "ci.*non-missing logical scalar"
  )
  expect_error(
    get_level_ts(res_tempssm, ci = 1),
    "ci.*non-missing logical scalar"
  )
})


test_that("ci_level is ignored when confidence intervals are not requested", {
  expect_no_error(
    get_level_ts(res_tempssm, ci = FALSE, ci_level = NA_real_)
  )
})


test_that("exo state information resolves names and positions", {
  fixture <- make_exo_accessor_fixture()

  info <- .tempssm_exo_state_info(fixture$res)

  expect_named(info, c("variables", "kfs", "alpha_hat", "indices"))
  expect_identical(info$variables, c("x1", "x2"))
  expect_identical(info$indices, c(1L, 2L))
  expect_identical(info$kfs, fixture$res$kfs)
  expect_identical(info$alpha_hat, fixture$res$kfs$alphahat)
})


test_that("exo state information returns NULL when coefficients are absent", {
  fixture <- make_exo_accessor_fixture()
  non_converged <- fixture$res
  non_converged$converged <- FALSE
  no_exogenous <- fixture$res
  no_exogenous$state_map$exogenous <- NULL

  expect_null(.tempssm_exo_state_info(non_converged))
  expect_null(.tempssm_exo_state_info(no_exogenous))
})


test_that("exo coefficient result construction preserves row correspondence", {
  fixture <- make_exo_accessor_fixture()
  state_info <- .tempssm_exo_state_info(fixture$res)
  ci_mat <- rbind(
    c(lwr = 0.1, upr = 0.4),
    c(lwr = -0.8, upr = -0.2)
  )

  result <- .new_exo_coef_result(state_info, ci_mat)

  expect_identical(
    result,
    data.frame(
      Variable = c("x1", "x2"),
      Coefficient = c(0.25, -0.5),
      lwr = c(0.1, -0.8),
      upr = c(0.4, -0.2)
    )
  )
})
