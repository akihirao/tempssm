test_that("tempssm rejects invalid model inputs before fitting", {
  expect_error(
    tempssm(1:10),
    "must be a <ts> object"
  )

  temp_ts <- ts(rnorm(12), frequency = 4)
  exo_matrix <- matrix(
    rnorm(12),
    ncol = 1,
    dimnames = list(NULL, "x")
  )

  expect_error(
    tempssm(temp_ts, exo_data = exo_matrix),
    "exo_data.*<ts> object"
  )
})


test_that("tempssm returns a well-formed failure object for fitting errors", {
  temp_ts <- ts(rnorm(12), frequency = 4)

  testthat::local_mocked_bindings(
    .define_build_model = function(...) {
      stop("synthetic fitting failure")
    },
    .package = "tempssm"
  )

  expect_warning(
    res <- tempssm(temp_ts),
    "Model fitting failed: synthetic fitting failure"
  )

  expect_s3_class(res, "tempssm")
  expect_false(res$converged)
  expect_null(res$model)
  expect_null(res$fit)
  expect_null(res$kfs)
  expect_identical(res$temp_data, temp_ts)
  expect_identical(res$call, quote(tempssm(temp_data = temp_ts)))
})


test_that("tempssm distinguishes non-convergence from fitting errors", {
  temp_ts <- ts(rnorm(12), frequency = 4)
  fit_count <- 0L

  testthat::local_mocked_bindings(
    fitSSM = function(model, inits, ...) {
      fit_count <<- fit_count + 1L
      list(
        model = model,
        optim.out = list(
          par = inits,
          convergence = if (fit_count == 2L) 1 else 0
        )
      )
    },
    KFS = function(model, ...) {
      n_states <- dim(model$T)[1]
      list(
        alphahat = matrix(
          0,
          nrow = length(temp_ts),
          ncol = n_states
        )
      )
    },
    .package = "tempssm"
  )

  expect_warning(
    res <- tempssm(temp_ts),
    "finished but did not converge"
  )

  expect_s3_class(res, "tempssm")
  expect_false(res$converged)
  expect_false(is.null(res$model))
  expect_false(is.null(res$fit))
  expect_false(is.null(res$kfs))
  expect_identical(res$call, quote(tempssm(temp_data = temp_ts)))
})
