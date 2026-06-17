# test-get_param_index.R

test_that(".get_param_index works with seasonal model", {

  idx <- .get_param_index(ar_order = 1, use_season = TRUE)

  expect_equal(idx$ar, 3)
  expect_equal(idx$var, 4)
  expect_equal(idx$H, 5)
})


test_that(".get_param_index works without seasonal model", {

  idx <- .get_param_index(ar_order = 1, use_season = FALSE)

  expect_equal(idx$ar, 2)
  expect_equal(idx$var, 3)
  expect_equal(idx$H, 4)
})


test_that("handles higher AR order correctly", {

  idx <- .get_param_index(ar_order = 3, use_season = TRUE)

  expect_equal(idx$ar, 3:5)
  expect_equal(idx$var, 6)
  expect_equal(idx$H, 7)
})


test_that("handles higher AR order without season", {

  idx <- .get_param_index(ar_order = 3, use_season = FALSE)

  expect_equal(idx$ar, 2:4)
  expect_equal(idx$var, 5)
  expect_equal(idx$H, 6)
})


test_that("returns correct structure", {

  idx <- .get_param_index(2, TRUE)

  expect_named(idx, c("ar", "var", "H"))
  expect_type(idx$ar, "integer")
  expect_length(idx$var, 1)
  expect_length(idx$H, 1)
})


test_that("ar_order = 1 returns scalar ar index", {

  idx <- .get_param_index(1, TRUE)

  expect_equal(length(idx$ar), 1)
})
