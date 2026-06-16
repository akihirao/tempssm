# tests/testthat/test-build_jma_url.R

test_that(".build_jma_url constructs correct URL", {
  url <- .build_jma_url("138")

  expect_equal(
    url,
    "https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/txt/area138.txt"
  )
})


test_that("numeric sea_area_id is handled", {
  url <- .build_jma_url(138)

  expect_equal(
    url,
    "https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/txt/area138.txt"
  )
})


test_that("character sea_area_id is preserved", {
  url <- .build_jma_url("001")

  expect_true(grepl("area001", url))
})


test_that("vector input returns vectorized output", {
  urls <- .build_jma_url(c("138", "139"))

  expect_length(urls, 2)
  expect_true(all(grepl("area", urls)))
})


test_that("handles empty input", {
  url <- .build_jma_url("")

  expect_true(grepl("area.txt", url))
})
