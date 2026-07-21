test_that("compare_ts_cv summarizes collected CV tables", {
  first <- tibble::tibble(
    fold = c(1, 2),
    converged = c(1, 0),
    MAE = c(1, NA),
    MASE_naive = c(0.5, NA),
    MASE_seasonal = c(0.7, NA)
  )
  second <- tibble::tibble(
    fold = c(1, 2),
    converged = c(1, 1),
    MAE = c(2, 4),
    MASE_naive = c(1, 3),
    MASE_seasonal = c(1.5, 2.5)
  )

  out <- compare_ts_cv(list(AR1 = first, AR2 = second))

  expect_s3_class(out, "tbl_df")
  expect_named(
    out,
    c(
      "model",
      "n_folds",
      "converged_n",
      "converged_rate",
      "mean_MAE",
      "mean_MASE_naive",
      "mean_MASE_seasonal"
    )
  )
  expect_identical(out$model, c("AR1", "AR2"))
  expect_identical(out$n_folds, c(2L, 2L))
  expect_identical(out$converged_n, c(1L, 2L))
  expect_identical(out$converged_rate, c(0.5, 1))
  expect_identical(out$mean_MAE, c(1, 3))
})


test_that("compare_ts_cv validates named list input", {
  tab <- tibble::tibble(
    fold = 1,
    converged = 1,
    MAE = 1,
    MASE_naive = 1,
    MASE_seasonal = 1
  )

  expect_error(
    compare_ts_cv(list(tab, tab)),
    "non-empty name"
  )

  expect_error(
    compare_ts_cv(list(AR1 = tab)),
    "at least two"
  )
})


test_that("compare_ts_cv validates table structure and folds", {
  tab <- tibble::tibble(
    fold = c(1, 2),
    converged = c(1, 1),
    MAE = c(1, 2),
    MASE_naive = c(1, 2),
    MASE_seasonal = c(1, 2)
  )
  missing_col <- tab[, setdiff(names(tab), "MASE_seasonal")]
  changed_folds <- tab
  changed_folds$fold <- c(1, 3)

  expect_error(
    compare_ts_cv(list(AR1 = tab, AR2 = missing_col)),
    "missing required column"
  )

  expect_error(
    compare_ts_cv(list(AR1 = tab, AR2 = changed_folds)),
    "same fold IDs"
  )
})
