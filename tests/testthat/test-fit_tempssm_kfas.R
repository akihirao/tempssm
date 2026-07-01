test_that(".fit_tempssm_kfas preserves the two-stage KFAS workflow", {
  calls <- new.env(parent = emptyenv())
  calls$fits <- list()
  calls$kfs <- NULL
  build_ssm <- list(id = "initial-model")
  updatefn <- function(pars, model) model
  inits <- c(-13, 0.5, -0.3, -5)
  withr::local_envvar(TEMPSSM_VERBOSITY = "none")

  testthat::local_mocked_bindings(
    fitSSM = function(model,
                      inits,
                      updatefn,
                      marginal,
                      method,
                      control) {
      stage <- length(calls$fits) + 1L
      calls$fits[[stage]] <- list(
        model = model,
        inits = inits,
        updatefn = updatefn,
        marginal = marginal,
        method = method,
        control = control
      )
      list(
        model = list(stage = stage),
        optim.out = list(
          par = inits + stage,
          convergence = 0
        )
      )
    },
    KFS = function(model, filtering, smoothing) {
      calls$kfs <- list(
        model = model,
        filtering = filtering,
        smoothing = smoothing
      )
      list(id = "kfs-result")
    },
    .package = "tempssm"
  )

  fitted <- .fit_tempssm_kfas(
    build_ssm = build_ssm,
    updatefn = updatefn,
    inits = inits,
    marginal = TRUE,
    maxit = 1000,
    reltol = 1e-10
  )

  expect_length(calls$fits, 2)
  expect_identical(calls$fits[[1]]$method, "Nelder-Mead")
  expect_identical(calls$fits[[2]]$method, "BFGS")
  expect_identical(calls$fits[[1]]$inits, inits)
  expect_identical(calls$fits[[2]]$inits, inits + 1)
  expect_true(calls$fits[[1]]$marginal)
  expect_true(calls$fits[[2]]$marginal)
  expect_identical(
    calls$fits[[1]]$control,
    list(maxit = 1000, reltol = 1e-10)
  )
  expect_identical(calls$fits[[2]]$control, calls$fits[[1]]$control)
  expect_identical(calls$kfs$model, fitted$fit$model)
  expect_identical(calls$kfs$filtering, c("state", "mean"))
  expect_identical(
    calls$kfs$smoothing,
    c("state", "mean", "disturbance")
  )
  expect_identical(fitted$kfs, list(id = "kfs-result"))
})


test_that(".fit_tempssm_kfas leaves fitting errors to its caller", {
  testthat::local_mocked_bindings(
    fitSSM = function(...) stop("synthetic KFAS failure"),
    .package = "tempssm"
  )

  expect_error(
    .fit_tempssm_kfas(
      build_ssm = list(),
      updatefn = function(pars, model) model,
      inits = c(-13, 0.5, -0.3, -5),
      marginal = FALSE,
      maxit = 1000,
      reltol = 1e-10
    ),
    "synthetic KFAS failure"
  )
})
