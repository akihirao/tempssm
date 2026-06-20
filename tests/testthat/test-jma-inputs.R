test_that("JMA SST retrieval functions validate sea area ID type", {
  expect_error(
    get_jma_sst_zoo(list(138)),
    "sea_area_id.*character or numeric"
  )

  expect_error(
    get_jma_sst_ts(list(138)),
    "sea_area_id.*character or numeric"
  )

  expect_error(
    get_jma_sst_ts(138, na_prop_max = "1"),
    "na_prop_max.*numeric"
  )
})
