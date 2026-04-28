#tests/testthat/test-diagnose_residuals.R

# JB test is not performed
test_that("diagnose_residuals returns a tibble", {

  diag <- diagnose_residuals(res_tempssm)

  expect_s3_class(diag, "tbl_df")
  expect_equal(nrow(diag), 1)

  expect_true(
    all(c("lb_stat", "lb_df", "lb_pvalue", "kurtosis") %in% colnames(diag))
  )
})


# JB test is performed
test_that("diagnose_residuals includes Jarque-Bera results when requested", {

  diag <- diagnose_residuals(res_tempssm, JB_test = TRUE)

  expect_s3_class(diag, "tbl_df")
  expect_true(
    all(c("jb_stat", "jb_pvalue") %in% colnames(diag))
  )
})


test_that("diagnose_residuals checks input class", {

  expect_error(
    diagnose_residuals(NULL),
    "`res` must be an object of class 'tempssm'"
  )
})
