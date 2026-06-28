make_jma_raw_fixture <- function(include_missing = FALSE) {
  jan_second <- if (include_missing) "-999" else "3"
  charToRaw(paste0(
    "Year,Month,Day,areaNo,flag,Temp\n",
    "2001,1,1,138,0,1\n",
    "2001,1,2,138,0,", jan_second, "\n",
    "2001,2,1,138,0,4\n",
    "2001,2,2,138,0,8\n",
    "2001,12,31,138,0,999\n",
    "END"
  ))
}


test_that("JMA retrieval functions share exact daily and monthly values", {
  testthat::local_mocked_bindings(
    .fetch_jma_raw = function(...) make_jma_raw_fixture(),
    .package = "tempssm"
  )

  daily <- get_jma_sst_zoo(138)
  monthly <- get_jma_sst_ts(138)

  expect_s3_class(daily, "zoo")
  expect_identical(as.numeric(zoo::coredata(daily)), c(1, 3, 4, 8))
  expect_identical(
    zoo::index(daily),
    as.Date(c("2001-01-01", "2001-01-02", "2001-02-01", "2001-02-02"))
  )
  expect_s3_class(monthly, "ts")
  expect_identical(as.numeric(monthly), c(2, 6))
  expect_identical(stats::start(monthly), c(2001, 1))
  expect_identical(stats::frequency(monthly), 12)
  expect_identical(colnames(monthly), "Temp")
})


test_that("get_jma_sst_ts passes na_prop_max to monthly aggregation", {
  testthat::local_mocked_bindings(
    .fetch_jma_raw = function(...) {
      make_jma_raw_fixture(include_missing = TRUE)
    },
    .package = "tempssm"
  )

  expect_warning(
    monthly <- get_jma_sst_ts("138", na_prop_max = 0.4),
    "More than 30% of aggregated monthly values are missing"
  )

  expect_true(is.na(as.numeric(monthly)[1]))
  expect_identical(as.numeric(monthly)[2], 6)
})


test_that("shared JMA retrieval builds the daily zoo representation", {
  testthat::local_mocked_bindings(
    .fetch_jma_raw = function(...) make_jma_raw_fixture(),
    .package = "tempssm"
  )

  result <- .retrieve_jma_sst_zoo("138")

  expect_s3_class(result, "zoo")
  expect_identical(colnames(result), "Temp")
  expect_identical(as.numeric(zoo::coredata(result)), c(1, 3, 4, 8))
})


test_that("shared JMA retrieval checks the daily missing-value ratio", {
  sparse_raw <- charToRaw(paste0(
    "Year,Month,Day,areaNo,flag,Temp\n",
    "2001,1,1,138,0,-999\n",
    "2001,1,2,138,0,-999\n",
    "2001,1,3,138,0,4\n",
    "2001,1,4,138,0,8\n",
    "2001,12,31,138,0,999\n",
    "END"
  ))
  testthat::local_mocked_bindings(
    .fetch_jma_raw = function(...) sparse_raw,
    .package = "tempssm"
  )

  expect_warning(
    result <- .retrieve_jma_sst_zoo("138"),
    "More than 30% of SST values are missing"
  )

  expect_identical(NROW(result), 4L)
})
