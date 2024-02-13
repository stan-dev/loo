library(loo)

context("pointwise convenience function")

loo1 <- suppressWarnings(loo(example_loglik_matrix()))

test_that("pointwise throws the right errors", {
  expect_error(
    pointwise(loo1, "xxx"),
    "'xxx' not found",
    fixed = TRUE
  )
  expect_error(
    pointwise(loo1, c("elpd_loo", "p_loo")),
    "length(estimate) == 1 is not TRUE",
    fixed = TRUE
  )
  expect_error(
    pointwise(loo1, 1),
    "is.character(estimate) is not TRUE",
    fixed = TRUE
  )
  loo1$pointwise <- NULL
  expect_error(
    pointwise(loo1, "xxx"),
    "No pointwise estimates found"
  )
})

test_that("pointwise returns correct estimate", {
  expect_equal(pointwise(loo1, "elpd_loo"), loo1$pointwise[, "elpd_loo"])
  expect_equal(pointwise(loo1, "mcse_elpd_loo"), loo1$pointwise[, "mcse_elpd_loo"])
  expect_equal(pointwise(loo1, "p_loo"), loo1$pointwise[, "p_loo"])
  expect_equal(pointwise(loo1, "looic"), loo1$pointwise[, "looic"])
  expect_equal(pointwise(loo1, "influence_pareto_k"), loo1$pointwise[, "influence_pareto_k"])
})
