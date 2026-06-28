test_that("JMA sea area validation preserves supported scalar IDs", {
  expect_identical(.prepare_jma_sea_area_id(138), 138)
  expect_identical(.prepare_jma_sea_area_id("001"), "001")
})


test_that("JMA sea area validation rejects unsupported inputs", {
  expect_error(
    .prepare_jma_sea_area_id(c(138, 139)),
    "sea_area_id.*length one"
  )
  expect_error(
    .prepare_jma_sea_area_id(list(138)),
    "sea_area_id.*character or numeric"
  )
})


test_that("JMA missing-value threshold validation is bounded", {
  expect_identical(.prepare_jma_na_prop_max(0), 0)
  expect_identical(.prepare_jma_na_prop_max(1), 1)
  expect_error(
    .prepare_jma_na_prop_max(NA_real_),
    "na_prop_max.*between 0 and 1"
  )
  expect_error(
    .prepare_jma_na_prop_max(Inf),
    "na_prop_max.*between 0 and 1"
  )
})
