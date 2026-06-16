# tests/testthat/test-compute_monthly_climatology.R

test_that("compute_monthly_climatology works for valid input", {
  temp_ts <- ts(
    rnorm(12 * 10, mean = 10),
    start = c(2000, 1),
    frequency = 12
  )

  res <- compute_monthly_climatology(temp_ts)

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 12)

  expect_named(res, c("Month", "Temperature"))
})


test_that("months are correctly ordered from 1 to 12", {
  temp_ts <- ts(
    rep(1:12, 10), # 明確な周期
    start = c(2000, 1),
    frequency = 12
  )

  res <- compute_monthly_climatology(temp_ts)

  expect_equal(res$Month, 1:12)
  expect_equal(res$Temperature, 1:12)
})


test_that("NA values are ignored in mean calculation", {
  temp_vec <- rep(1:12, 5)
  temp_vec[1] <- NA # Januaryの1つをNAに

  temp_ts <- ts(temp_vec, frequency = 12)

  res <- compute_monthly_climatology(temp_ts)

  expect_false(is.na(res$Temperature[1]))
})


test_that("returns NA when all values for a month are NA", {
  temp_vec <- rep(1:12, 5)
  temp_vec[seq(1, length(temp_vec), by = 12)] <- NA # JanuaryすべてNA

  temp_ts <- ts(temp_vec, frequency = 12)

  res <- compute_monthly_climatology(temp_ts)

  expect_true(is.na(res$Temperature[1]))
})


test_that("errors when input is not ts", {
  expect_error(
    compute_monthly_climatology(rnorm(10)),
    "must be a"
  )
})


test_that("errors when frequency is not 12", {
  temp_ts <- ts(rnorm(10), frequency = 4)

  expect_error(
    compute_monthly_climatology(temp_ts),
    "must be monthly"
  )
})
