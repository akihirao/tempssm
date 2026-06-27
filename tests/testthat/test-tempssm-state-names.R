test_that(".make_tempssm_state_names handles all model configurations", {
  expect_identical(
    .make_tempssm_state_names(NULL, 4, FALSE, 1, 3),
    c("level", "slope", "arima1")
  )

  expect_identical(
    .make_tempssm_state_names(NULL, 4, TRUE, 1, 6),
    c(
      "level",
      "slope",
      "sea_dummy1",
      "sea_dummy2",
      "sea_dummy3",
      "arima1"
    )
  )

  expect_identical(
    .make_tempssm_state_names(c("pdo", "enso"), 4, FALSE, 1, 5),
    c("pdo", "enso", "level", "slope", "arima1")
  )

  expect_identical(
    .make_tempssm_state_names(c("pdo", "enso"), 4, TRUE, 2, 9),
    c(
      "pdo",
      "enso",
      "level",
      "slope",
      "sea_dummy1",
      "sea_dummy2",
      "sea_dummy3",
      "arima1",
      "arima2"
    )
  )
})


test_that(".make_tempssm_state_names detects KFAS state-count mismatches", {
  expect_error(
    .make_tempssm_state_names(NULL, 4, TRUE, 1, 5),
    "produced 6 names.*state matrix has 5 columns"
  )
})
