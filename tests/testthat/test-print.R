# test-print.R

test_that("print.tempssm works without error", {
  expect_no_error(print(res_tempssm))
})


test_that("print.tempssm returns input object invisibly", {

  out <- print(res_tempssm)
  expect_identical(out, res_tempssm)
})


test_that("print.tempssm handles NULL exogenous_data", {

  res_tempssm$exogenous_data <- NULL

  expect_no_error(print(res_tempssm))
})


test_that("print.tempssm outputs expected text", {


  output <- capture.output(print(res_tempssm))

  expect_true(any(grepl("tempssm model fit", output)))
})



test_that("print.tempssm handles non-converged case", {

  res <- res_tempssm
  res$fit$optim.out$convergence <- 1  # 強制的に非収束

  output <- capture.output(print(res))

  expect_true(any(grepl("Converged", output)))
})



test_that("print.tempssm displays time series info", {

  output <- capture.output(print(res_tempssm))

  expect_true(any(grepl("Length", output)))
  expect_true(any(grepl("Frequency", output)))
  expect_true(any(grepl("Start / End", output)))
})



test_that("print.tempssm includes logLik output", {

  output <- capture.output(print(res_tempssm))

  expect_true(any(grepl("LogLik", output)))
})



test_that("print.tempssm explicitly shows NULL exogenous", {

  res <- res_tempssm
  res$exogenous_data <- NULL

  output <- capture.output(print(res))

  expect_true(any(grepl("NULL", output)))
})





