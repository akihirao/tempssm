# test-print.R

test_that("print.tempssm basic output structure", {

  output <- capture.output(print(res_tempssm))

  expect_true(any(grepl("tempssm model fit", output)))
  expect_true(any(grepl("Data:", output)))
  expect_true(any(grepl("Optimization:", output)))
  expect_true(any(grepl("Use summary", output)))
})


test_that("print.tempssm displays correct time series info", {

  output <- capture.output(print(res_tempssm))

  expect_true(any(grepl(paste(length(res_tempssm$temp_data)), output)))
  expect_true(any(grepl(paste(frequency(res_tempssm$temp_data)), output)))
})


test_that("print.tempssm handles exogenous data correctly", {

  res <- res_tempssm_exo  # setupで用意

  output <- capture.output(print(res))

  expect_true(any(grepl("Exogenous variable", output)))
  expect_true(any(grepl("No. variables", output)))
})


test_that("print.tempssm shows NULL exogenous correctly", {

  res <- res_tempssm
  res$exogenous_data <- NULL

  output <- capture.output(print(res))

  expect_true(any(grepl("Exogenous variables: NULL", output)))
})


test_that("print.tempssm reflects convergence status", {

  res <- res_tempssm
  res$fit$optim.out$convergence <- 1

  output <- capture.output(print(res))

  expect_true(any(grepl("FALSE", output)))
})


test_that("print.tempssm prints logLik as numeric", {

  output <- capture.output(print(res_tempssm))

  expect_true(any(grepl("LogLik", output)))

  # 数値っぽい形式チェック
  expect_true(any(grepl("LogLik.*[0-9]", output)))
})


test_that("print.tempssm does not modify object", {

  res_copy <- res_tempssm

  print(res_copy)

  expect_identical(res_copy, res_tempssm)
})


test_that("print.tempssm returns invisibly", {

  expect_invisible(out <- print(res_tempssm))
  expect_identical(out, res_tempssm)
})
