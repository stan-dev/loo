library(loo)

context("pointwise convenience function")

loo1 <- suppressWarnings(loo(example_loglik_matrix()))

test_that("pointwise throws the right errors", {
  expect_error(
    pointwise(loo1, "xxx"),
    "'xxx' not found"
  )
  loo1$pointwise <- NULL
  expect_error(
    pointwise(loo1, "xxx"),
    "No pointwise estimates found"
  )
})

test_that("pointwise returns matrix if no estimate specified", {
  pw <- pointwise(loo1)
  expect_true(is.matrix(pw))
  expect_equal(pw, loo1$pointwise)
})

test_that("pointwise returns correct single estimate", {
  expect_equal(pointwise(loo1, "elpd_loo"), loo1$pointwise[, "elpd_loo"])
  expect_equal(pointwise(loo1, "mcse_elpd_loo"), loo1$pointwise[, "mcse_elpd_loo"])
  expect_equal(pointwise(loo1, "p_loo"), loo1$pointwise[, "p_loo"])
  expect_equal(pointwise(loo1, "looic"), loo1$pointwise[, "looic"])
  expect_equal(pointwise(loo1, "influence_pareto_k"), loo1$pointwise[, "influence_pareto_k"])
})
