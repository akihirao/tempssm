# test-tempssm-inits.R

test_that(".tempssm_inits_length matches model structure", {
  expect_identical(.tempssm_inits_length(1, TRUE), 5L)
  expect_identical(.tempssm_inits_length(1, FALSE), 4L)
  expect_identical(.tempssm_inits_length(3, TRUE), 7L)
  expect_identical(.tempssm_inits_length(3, FALSE), 6L)
})


test_that(".default_tempssm_inits follows the parameter order", {
  expect_identical(
    .default_tempssm_inits(ar_order = 1, use_season = TRUE),
    c(-13, -7, 0.5, -0.3, -5)
  )

  expect_identical(
    .default_tempssm_inits(ar_order = 3, use_season = FALSE),
    c(-13, 0.5, 0, 0, -0.3, -5)
  )
})


test_that(".prepare_tempssm_inits returns valid supplied values unchanged", {
  supplied <- c(-10, -6, 0.2, -0.1, -4)

  expect_identical(
    .prepare_tempssm_inits(supplied, ar_order = 1, use_season = TRUE),
    supplied
  )
})


test_that(".prepare_tempssm_inits supplies and validates defaults", {
  prepared <- .prepare_tempssm_inits(NULL, ar_order = 2, use_season = FALSE)

  expect_identical(prepared, c(-13, 0.5, 0, -0.3, -5))
})


test_that(".prepare_tempssm_inits rejects invalid supplied values", {
  expect_error(
    .prepare_tempssm_inits(c(1, 2), ar_order = 1, use_season = TRUE),
    "inits.*length"
  )

  expect_error(
    .prepare_tempssm_inits(rep("x", 5), ar_order = 1, use_season = TRUE),
    "inits.*numeric"
  )
})
