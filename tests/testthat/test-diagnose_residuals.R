# tests/testthat/test-diagnose_residuals.R

# JB test is not performed
test_that("diagnose_residuals returns a tibble", {
  diag <- diagnose_residuals(res_tempssm)

  expect_s3_class(diag, "tbl_df")
  expect_identical(nrow(diag), 1L)

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
    "`res` must be an object of class <tempssm>."
  )
})


test_that("diagnose_residuals validates scalar argument lengths", {
  expect_error(
    diagnose_residuals(res_tempssm, JB_test = c(TRUE, FALSE)),
    "JB_test.*length one"
  )
})


test_that("diagnose_residuals validates scalar argument types", {
  expect_error(
    diagnose_residuals(res_tempssm, JB_test = 1),
    "JB_test.*logical"
  )
})


test_that(".kurtosis removes missing values when requested", {
  x <- c(1, 2, NA, 3, 4)

  expect_true(is.na(.kurtosis(x)))
  expect_false(is.na(.kurtosis(x, na.rm = TRUE)))
})


test_that("plot_tempssm_residual_diagnostics can save plots", {
  prefix <- file.path(tempdir(), "tempssm-residuals")
  check_file <- paste0(prefix, "_check.png")
  qq_file <- paste0(prefix, "_qq.png")

  withr::defer(unlink(c(check_file, qq_file)))

  expect_invisible(
    plot_tempssm_residual_diagnostics(res_tempssm, save = TRUE, prefix = prefix)
  )
  expect_true(file.exists(check_file))
  expect_true(file.exists(qq_file))
})


test_that("plot_tempssm_residual_diagnostics normalizes prefix extensions", {
  prefix <- file.path(tempdir(), "tempssm-residuals.png")
  check_file <- file.path(tempdir(), "tempssm-residuals_check.png")
  qq_file <- file.path(tempdir(), "tempssm-residuals_qq.png")
  unexpected_file <- paste0(prefix, "_check.png")

  withr::defer(unlink(c(check_file, qq_file, unexpected_file)))

  expect_invisible(
    plot_tempssm_residual_diagnostics(res_tempssm, save = TRUE, prefix = prefix)
  )
  expect_true(file.exists(check_file))
  expect_true(file.exists(qq_file))
  expect_false(file.exists(unexpected_file))
})


test_that(
  "plot_tempssm_residual_diagnostics validates scalar argument lengths",
  {
  expect_error(
    plot_tempssm_residual_diagnostics(res_tempssm, save = c(TRUE, FALSE)),
    "save.*length one"
  )

  expect_error(
    plot_tempssm_residual_diagnostics(
      res_tempssm,
      prefix = c("a", "b")
    ),
    "prefix.*length one"
  )
  }
)


test_that("plot_tempssm_residual_diagnostics validates scalar argument types", {
  expect_error(
    plot_tempssm_residual_diagnostics(res_tempssm, save = 1),
    "save.*logical"
  )

  expect_error(
    plot_tempssm_residual_diagnostics(res_tempssm, prefix = 1),
    "prefix.*character"
  )
})
