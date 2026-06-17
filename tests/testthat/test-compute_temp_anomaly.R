# test-compute_temp_anomaly.R


test_that("compute_temp_anomaly works without baseline", {
  temp_ts <- ts(
    rnorm(12 * 10, mean = 10),
    start = c(2000, 1),
    frequency = 12
  )

  anom <- compute_temp_anomaly(temp_ts)

  expect_s3_class(anom, "ts")
  expect_equal(frequency(anom), 12)
  expect_equal(start(anom), start(temp_ts))
  expect_length(anom, length(temp_ts))
})


test_that("compute_temp_anomaly works with baseline", {
  temp_ts <- ts(
    rnorm(12 * 20, mean = 10),
    start = c(2000, 1),
    frequency = 12
  )

  anom <- compute_temp_anomaly(temp_ts, baseline = c(2005, 2010))

  expect_s3_class(anom, "ts")
  expect_length(anom, length(temp_ts))
})


test_that("anomalies have near-zero monthly mean", {
  set.seed(123)

  temp_ts <- ts(
    rnorm(12 * 30, mean = 15),
    start = c(1990, 1),
    frequency = 12
  )

  anom <- compute_temp_anomaly(temp_ts)

  monthly_means <- tapply(
    as.numeric(anom),
    cycle(anom),
    mean
  )

  expect_true(all(abs(monthly_means) < 1e-10))
})


test_that("non-ts input triggers error", {
  expect_error(
    compute_temp_anomaly(rnorm(10)),
    "must be a"
  )
})


test_that("non-monthly ts triggers error", {
  temp_ts <- ts(rnorm(10), frequency = 4)

  expect_error(
    compute_temp_anomaly(temp_ts),
    "monthly"
  )
})


test_that("baseline with no data triggers error", {
  temp_ts <- ts(
    rnorm(24),
    start = c(2000, 1),
    frequency = 12
  )

  expect_error(
    compute_temp_anomaly(temp_ts, baseline = c(1980, 1985)),
    "No data available"
  )
})


test_that("output values are numeric", {
  temp_ts <- ts(rnorm(24), frequency = 12)

  anom <- compute_temp_anomaly(temp_ts)

  expect_true(is.numeric(as.numeric(anom)))
})
