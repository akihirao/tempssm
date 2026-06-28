# tests/testthat/test-daily_zoo_to_monthly_ts.R

test_that("daily_zoo_to_monthly_ts works for valid input", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 60)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(60)),
    order.by = dates
  )

  res <- daily_zoo_to_monthly_ts(zoo_obj)

  expect_s3_class(res, "ts")
  expect_identical(frequency(res), 12)
  expect_gt(length(res), 2) # 2 months
  expect_identical(colnames(res), "Temp")
})


test_that("monthly means and the missing-value threshold are exact", {
  dates <- as.Date(c(
    "2001-01-01", "2001-01-02", "2001-01-03", "2001-01-04",
    "2001-02-01", "2001-02-02", "2001-02-03", "2001-02-04"
  ))
  zoo_obj <- zoo::zoo(
    data.frame(Temp = c(1, 3, NA, 5, 2, 4, 6, 8)),
    order.by = dates
  )

  at_threshold <- daily_zoo_to_monthly_ts(
    zoo_obj,
    na_prop_max = 0.25
  )
  expect_warning(
    below_threshold <- daily_zoo_to_monthly_ts(
      zoo_obj,
      na_prop_max = 0.24
    ),
    "More than 30%"
  )

  expect_identical(as.numeric(at_threshold), c(3, 5))
  expect_true(is.na(as.numeric(below_threshold)[1]))
  expect_identical(as.numeric(below_threshold)[2], 5)
})


test_that("na.rm FALSE preserves missing monthly means", {
  dates <- as.Date(c(
    "2001-01-01", "2001-01-02",
    "2001-02-01", "2001-02-02"
  ))
  zoo_obj <- zoo::zoo(
    data.frame(Temp = c(1, NA, 3, 5)),
    order.by = dates
  )

  expect_warning(
    result <- daily_zoo_to_monthly_ts(zoo_obj, na.rm = FALSE),
    "More than 30%"
  )

  expect_true(is.na(as.numeric(result)[1]))
  expect_identical(as.numeric(result)[2], 4)
})


test_that("POSIXct indices are aggregated by calendar month", {
  dates <- as.POSIXct(
    c("2001-01-01", "2001-01-02", "2001-02-01", "2001-02-02"),
    tz = "UTC"
  )
  zoo_obj <- zoo::zoo(
    data.frame(Temp = c(1, 3, 4, 8)),
    order.by = dates
  )

  result <- daily_zoo_to_monthly_ts(zoo_obj)

  expect_identical(as.numeric(result), c(2, 6))
  expect_identical(stats::start(result), c(2001, 1))
  expect_identical(stats::frequency(result), 12)
})


test_that("var argument selects correct column", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(A = rnorm(30), B = rnorm(30)),
    order.by = dates
  )

  res <- daily_zoo_to_monthly_ts(zoo_obj, var = "B")

  expect_identical(colnames(res), "B")
})


test_that("errors when variable not found", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(30)),
    order.by = dates
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, var = "X"),
    "not found"
  )
})


test_that("errors for invalid na_prop_max", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(30)),
    order.by = dates
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, na_prop_max = -1),
    "na_prop_max.*between 0 and 1"
  )
  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, na_prop_max = 2),
    "na_prop_max.*between 0 and 1"
  )
  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, na_prop_max = c(0.2, 0.5)),
    "na_prop_max.*length one"
  )
})


test_that("daily_zoo_to_monthly_ts validates scalar argument lengths", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(30)),
    order.by = dates
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, var = c("Temp", "Other")),
    "var.*length one"
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, na.rm = c(TRUE, FALSE)),
    "na.rm.*length one"
  )
})


test_that("daily_zoo_to_monthly_ts validates scalar argument types", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(30)),
    order.by = dates
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, var = 1),
    "var.*character"
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, na.rm = 1),
    "na.rm.*logical"
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj, na_prop_max = "1"),
    "na_prop_max.*numeric"
  )
})


test_that("errors when index is not Date/POSIXt", {
  zoo_obj <- zoo::zoo(
    data.frame(Temp = rnorm(10)),
    order.by = 1:10
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj),
    "Index of the zoo object"
  )
})


test_that("errors when zoo index is not strictly increasing", {
  zoo_obj <- suppressWarnings(
    zoo::zoo(
      data.frame(Temp = rnorm(3)),
      order.by = as.Date(c("2001-01-01", "2001-01-01", "2001-01-02"))
    )
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj),
    "strictly increasing"
  )
})


test_that("na.rm works correctly in aggregation", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)

  vals <- c(1, 2, NA, 3, 4, NA, 5, rep(1, 23))

  zoo_obj <- zoo::zoo(data.frame(Temp = vals), order.by = dates)

  res <- daily_zoo_to_monthly_ts(zoo_obj, na.rm = TRUE)

  expect_false(is.na(res[1]))
})


test_that("undefined values after aggregation trigger error", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 30)
  zoo_obj <- zoo::zoo(
    data.frame(Temp = c(Inf, rep(1, 29))),
    order.by = dates
  )

  expect_error(
    daily_zoo_to_monthly_ts(zoo_obj),
    "Inf"
  )
})


test_that("months without daily observations are explicit NA values", {
  dates <- c(
    seq.Date(as.Date("2001-01-01"), by = "day", length.out = 5),
    seq.Date(as.Date("2001-03-01"), by = "day", length.out = 5),
    seq.Date(as.Date("2001-04-01"), by = "day", length.out = 5)
  )

  zoo_obj <- zoo::zoo(
    data.frame(Temp = seq_along(dates)),
    order.by = dates
  )

  expect_warning(
    res <- daily_zoo_to_monthly_ts(zoo_obj),
    "Implicit missing months"
  )

  expect_identical(start(res), c(2001, 1))
  expect_length(res, 4)
  expect_true(is.na(as.numeric(res)[2]))
})


test_that("warns when many NA in aggregated result", {
  dates <- seq.Date(as.Date("2001-01-01"), by = "day", length.out = 60)

  vals <- c(rep(NA, 40), rnorm(20)) # 66% NA

  zoo_obj <- zoo::zoo(data.frame(Temp = vals), order.by = dates)

  expect_warning(
    daily_zoo_to_monthly_ts(zoo_obj),
    "More than 30%"
  )
})
