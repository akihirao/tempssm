test_that("trim_ts_overlap preserves exact overlap values and missing data", {
  temp_data <- ts(seq_len(12), start = c(2000, 1), frequency = 4)
  temp_data[7] <- NA_real_
  exo_data <- ts(
    cbind(x = 101:108, y = 201:208),
    start = c(2001, 1),
    frequency = 4
  )

  result <- trim_ts_overlap(
    temp_data,
    exo_data,
    temp_name = "temperature",
    exo_name = c("x", "y")
  )

  expect_identical(as.numeric(result$temperature), temp_data[5:12])
  expect_true(is.na(as.numeric(result$temperature)[3]))
  expect_identical(
    as.numeric(result$exogenous[, "x"]),
    as.numeric(101:108)
  )
  expect_identical(
    as.numeric(result$exogenous[, "y"]),
    as.numeric(201:208)
  )
  expect_identical(stats::start(result$temperature), c(2001, 1))
  expect_identical(stats::end(result$temperature), c(2002, 4))
})


test_that("trim_ts_overlap supports matching non-monthly frequencies", {
  for (freq in c(4, 24)) {
    temp_data <- ts(
      seq_len(freq * 2),
      start = c(2000, 1),
      frequency = freq
    )
    exo_data <- ts(
      seq_len(freq),
      start = c(2001, 1),
      frequency = freq
    )

    result <- trim_ts_overlap(
      temp_data,
      exo_data,
      exo_name = "x"
    )

    expect_identical(stats::frequency(result$temperature), freq)
    expect_identical(stats::frequency(result$exogenous), freq)
    expect_identical(NROW(result$temperature), as.integer(freq))
  }
})


test_that("trim_ts_overlap rejects different frequencies", {
  temp_data <- ts(seq_len(8), start = c(2000, 1), frequency = 4)
  exo_data <- ts(seq_len(24), start = c(2000, 1), frequency = 12)

  expect_error(
    trim_ts_overlap(temp_data, exo_data, exo_name = "x"),
    "same frequency|not all series have the same frequency"
  )
})


test_that("trim_ts_overlap rejects non-overlapping series", {
  temp_data <- ts(seq_len(4), start = c(2000, 1), frequency = 4)
  exo_data <- ts(seq_len(4), start = c(2002, 1), frequency = 4)

  suppressWarnings(
    expect_error(
      trim_ts_overlap(temp_data, exo_data, exo_name = "x"),
      "No overlapping time period"
    )
  )
})


test_that("trim overlap alias resolution preserves current and legacy inputs", {
  current <- ts(seq_len(4), frequency = 4)
  legacy <- ts(5:8, frequency = 4)

  expect_identical(
    .resolve_trim_ts_alias(
      current,
      legacy = NULL,
      current_name = "current",
      legacy_name = "legacy"
    ),
    current
  )
  expect_identical(
    .resolve_trim_ts_alias(
      legacy = legacy,
      current_name = "current",
      legacy_name = "legacy"
    ),
    legacy
  )
  expect_error(
    .resolve_trim_ts_alias(
      current,
      legacy = legacy,
      current_name = "current",
      legacy_name = "legacy"
    ),
    "Use either.*current.*legacy"
  )
})


test_that("trim overlap name preparation validates and supplies names", {
  expect_warning(
    default_names <- .prepare_trim_ts_names("temp", NULL, 2L),
    "exo_name.*not supplied"
  )
  supplied_names <- .prepare_trim_ts_names(
    "temperature",
    c("x", "y"),
    2L
  )

  expect_identical(default_names$temp_name, "temp")
  expect_identical(default_names$exo_name, c("var1", "var2"))
  expect_identical(supplied_names$temp_name, "temperature")
  expect_identical(supplied_names$exo_name, c("x", "y"))
  expect_error(
    .prepare_trim_ts_names("temp", c("x", NA_character_), 2L),
    "exo_name.*must not contain missing values"
  )
})


test_that("trim overlap intersection returns separated exact series", {
  temp_data <- ts(seq_len(8), start = c(2000, 1), frequency = 4)
  exo_data <- ts(
    cbind(x = 11:14, y = 21:24),
    start = c(2001, 1),
    frequency = 4
  )

  overlap <- .intersect_trim_ts(temp_data, exo_data)

  expect_named(overlap, c("temperature", "exogenous"))
  expect_identical(as.numeric(overlap$temperature), as.numeric(5:8))
  expect_identical(
    as.numeric(overlap$exogenous[, 1L]),
    as.numeric(11:14)
  )
  expect_identical(
    as.numeric(overlap$exogenous[, 2L]),
    as.numeric(21:24)
  )
})


test_that("trim overlap result construction applies requested labels", {
  overlap <- list(
    temperature = ts(matrix(1:4, ncol = 1), frequency = 4),
    exogenous = ts(cbind(x = 11:14, y = 21:24), frequency = 4)
  )

  result <- .new_trimmed_ts_result(
    overlap,
    temp_name = "temperature",
    exo_name = c("wind", "pressure")
  )

  expect_named(result, c("temperature", "exogenous"))
  expect_identical(colnames(result$temperature), "temperature")
  expect_identical(colnames(result$exogenous), c("wind", "pressure"))
  expect_identical(stats::frequency(result$temperature), 4)
  expect_identical(stats::frequency(result$exogenous), 4)
})
