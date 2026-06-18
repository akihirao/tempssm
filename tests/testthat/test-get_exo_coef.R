# tests/testthat/test-get_exo_coef.R

test_that("get_exo_coef returns data.frame with correct structure", {
  if (is.null(res_tempssm_exo$exogenous_data)) skip("no exogenous variables")

  df <- get_exo_coef(res_tempssm_exo)

  expect_s3_class(df, "data.frame")

  expect_named(
    df,
    c("Variable", "Coefficient", "lwr", "upr")
  )
})


test_that("get_exo_coef returns numeric columns", {
  if (is.null(res_tempssm_exo$exogenous_data)) skip()

  df <- get_exo_coef(res_tempssm_exo)

  expect_type(df$Coefficient, "double")
  expect_type(df$lwr, "double")
  expect_type(df$upr, "double")
})


test_that("returns NULL when model did not converge", {
  bad_res <- res_tempssm_exo
  bad_res$converged <- FALSE

  expect_null(get_exo_coef(bad_res))
})


test_that("coefficient matches alphahat", {
  if (is.null(res_tempssm_exo$exogenous_data)) skip()

  df <- get_exo_coef(res_tempssm_exo)

  var_name <- df$Variable

  alpha <- res_tempssm_exo$kfs$alphahat

  expect_identical(
    df$Coefficient,
    as.numeric(alpha[1, var_name])
  )
})


test_that("get_exo_coef checks inputs", {
  expect_error(
    get_exo_coef(NULL),
    "`res` must be an object of class"
  )

  expect_error(
    get_exo_coef(res_tempssm_exo, ci_level = 2),
    "`ci_level` must be a numeric value between 0 and 1"
  )
})
