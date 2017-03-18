library(loo)

context("loo_expectation")

test_that("loo_expectation return types are correct", {
  set.seed(123)
  x <- matrix(rnorm(200), 20, 10)
  lw <- matrix(rnorm(200), 20, 10)

  expect_null(dim(loo_expectation(x[,1], lw[,1], type="mean")))
  expect_equal(length(loo_expectation(x[,1], lw[,1], type="var")), 1)
  expect_equal(
    length(loo_expectation(x[,1], lw[,1], type="quantile", probs = c(.25, .75))),
    2
  )

  expect_null(dim(loo_expectation(x, lw, type="mean")))
  expect_null(dim(loo_expectation(x, lw, type="var")))
  expect_equal(length(loo_expectation(x, lw, type="mean")), 10)
  expect_equal(length(loo_expectation(x, lw, type="var")), 10)
  expect_equal(
    dim(loo_expectation(x, lw, type="quantile", probs = c(.25, .75))),
    c(2, 10)
  )
})

test_that("loo_expectation.default equal to reference", {
  set.seed(124)
  x <- rnorm(100)
  lw <- rnorm(100)

  expect_equal_to_reference(
    loo_expectation(x, lw, type = "mean"),
    "loo_expectation_default_mean.rds"
  )
  expect_equal_to_reference(
    loo_expectation(x, lw, type = "var"),
    "loo_expectation_default_var.rds"
  )
  expect_equal_to_reference(
    loo_expectation(x, lw, type = "quantile", probs = 0.5),
    "loo_expectation_default_quantile_50.rds"
  )
  expect_equal_to_reference(
    loo_expectation(x, lw, type = "quantile", probs = c(0.1, 0.9)),
    "loo_expectation_default_quantile_10_90.rds"
  )
})

test_that("loo_expectation.matrix equal to reference", {
  set.seed(125)
  x <- matrix(rnorm(200), 20, 10)
  lw <- matrix(rnorm(200), 20, 10)

  expect_equal_to_reference(
    loo_expectation(x, lw, type = "mean"),
    "loo_expectation_matrix_mean.rds"
  )
  expect_equal_to_reference(
    loo_expectation(x, lw, type = "var"),
    "loo_expectation_matrix_var.rds"
  )
  expect_equal_to_reference(
    loo_expectation(x, lw, type = "quantile", probs = 0.5),
    "loo_expectation_matrix_quantile_50.rds"
  )
  expect_equal_to_reference(
    loo_expectation(x, lw, type = "quantile", probs = c(0.1, 0.9)),
    "loo_expectation_matrix_quantile_10_90.rds"
  )
})

test_that("loo_expectation throws correct errors", {
  expect_error(loo_expectation(1, 1, type = "purple"))
  expect_error(
    loo_expectation(1, 1, type = "quantile", probs = 2),
    "all(probs > 0 & probs < 1) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    loo_expectation("a", "b"),
    "is.numeric(x) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    loo_expectation(1:10, 1:5),
    "length(lw) == length(x) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    loo_expectation(cbind(1:10, 1:10), cbind(1:5, 1:5)),
    "identical(dim(x), dim(lw)) is not TRUE",
    fixed = TRUE
  )
})

