# tests/testthat/test-utils-verbosity.R

test_that(".get_verbosity prefers environment variable", {
  withr::local_envvar(TEMPSSM_VERBOSITY = "debug")
  withr::local_options(list(tempssm.verbosity = "none"))

  expect_identical(.get_verbosity(), "debug")
})


test_that(".get_verbosity falls back to option and validates values", {
  withr::local_envvar(TEMPSSM_VERBOSITY = "")
  withr::local_options(list(tempssm.verbosity = "none"))

  expect_identical(.get_verbosity(), "none")

  withr::local_options(list(tempssm.verbosity = "verbose"))

  expect_identical(.get_verbosity(), "inform")
})


test_that("verbosity helpers respect none and debug modes", {
  withr::local_envvar(TEMPSSM_VERBOSITY = "none")

  expect_silent(.tempssm_cli_inform("hidden message"))
  expect_silent(.tempssm_cli_debug("hidden debug message"))

  withr::local_envvar(TEMPSSM_VERBOSITY = "debug")

  expect_message(.tempssm_cli_debug("shown debug message"), "shown debug")
})
