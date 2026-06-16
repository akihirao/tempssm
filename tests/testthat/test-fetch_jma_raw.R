# tests/testthat/test-fetch_jma_raw.R

test_that(".fetch_jma_raw returns raw vector", {
  fake_raw <- charToRaw("dummy data")

  mock_resp <- list()

  mockery::stub(
    .fetch_jma_raw,
    "httr2::req_perform",
    function(...) mock_resp
  )

  mockery::stub(
    .fetch_jma_raw,
    "httr2::resp_body_raw",
    function(...) fake_raw
  )

  res <- .fetch_jma_raw("138")

  expect_type(res, "raw")
  expect_equal(res, fake_raw)
})


test_that("numeric sea_area_id is accepted", {
  fake_raw <- charToRaw("ok")

  mockery::stub(
    .fetch_jma_raw,
    "httr2::req_perform",
    function(...) list()
  )

  mockery::stub(
    .fetch_jma_raw,
    "httr2::resp_body_raw",
    function(...) fake_raw
  )

  expect_no_error(
    .fetch_jma_raw(138)
  )
})


test_that("fails on HTTP error", {
  mockery::stub(
    .fetch_jma_raw,
    "httr2::req_perform",
    function(...) stop("network error")
  )

  expect_error(
    .fetch_jma_raw("138"),
    "Failed to access JMA SST endpoint"
  )
})
