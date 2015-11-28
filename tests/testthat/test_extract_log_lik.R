library(loo)

context("extract_log_lik")
test_that("extract_log_lik throws appropriate errors", {
  x1 <- rnorm(100)
  expect_error(extract_log_lik(x1), regexp = "Not a stanfit object")
  x2 <- structure(x1, class = "stanfit")
  expect_error(extract_log_lik(x2), "not an S4 object")
})
