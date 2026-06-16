# tests/testthat/test-build_sst_zoo.R

test_that(".build_sst_zoo returns zoo object", {
  sst <- data.frame(
    date = as.Date(c("2000-01-01", "2000-01-02")),
    Temp = c(15.2, 16.0)
  )

  res <- .build_sst_zoo(sst)

  expect_s3_class(res, "zoo")
})


test_that("zoo object has correct column name", {
  sst <- data.frame(
    date = as.Date("2000-01-01") + 0:1,
    Temp = c(10, 11)
  )

  res <- .build_sst_zoo(sst)

  expect_named(res, "Temp")
})


test_that("index matches input dates", {
  dates <- as.Date(c("2000-01-01", "2000-01-02"))

  sst <- data.frame(
    date = dates,
    Temp = c(1, 2)
  )

  res <- .build_sst_zoo(sst)

  expect_equal(zoo::index(res), dates)
})


test_that("values are preserved", {
  temp <- c(15.5, 16.2)

  sst <- data.frame(
    date = as.Date("2000-01-01") + 0:1,
    Temp = temp
  )

  res <- .build_sst_zoo(sst)

  expect_equal(as.numeric(res), temp)
})


test_that("NA values are preserved", {
  sst <- data.frame(
    date = as.Date("2000-01-01") + 0:2,
    Temp = c(1, NA, 3)
  )

  res <- .build_sst_zoo(sst)

  expect_true(is.na(res[2]))
})


test_that("works with single observation", {
  sst <- data.frame(
    date = as.Date("2000-01-01"),
    Temp = 10
  )

  res <- .build_sst_zoo(sst)

  expect_equal(length(res), 1)
})


test_that("handles empty input", {
  sst <- data.frame(
    date = as.Date(character(0)),
    Temp = numeric(0)
  )

  res <- .build_sst_zoo(sst)

  expect_s3_class(res, "zoo")
  expect_equal(length(res), 0)
})
