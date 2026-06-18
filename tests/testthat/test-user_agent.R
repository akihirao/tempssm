# tests/testthat/test-user_agent.R

test_that(".user_agent returns default value when option not set", {
  withr::local_options(list(tempssm.user_agent = NULL))

  ua <- .user_agent()

  expect_type(ua, "character")
  expect_length(ua, 1)

  expect_match(ua, "^tempssm/")
  expect_match(ua, "github.com/yourname/tempssm")
})


test_that("default user agent contains package version", {
  withr::local_options(list(tempssm.user_agent = NULL))

  ua <- .user_agent()

  pkg_ver <- as.character(utils::packageVersion("tempssm"))

  expect_match(ua, pkg_ver)
})
