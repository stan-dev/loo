library(loo)
options(loo.cores = 1)

set.seed(123)
x <- matrix(rnorm(5000), 100, 50)

context("psislw")
test_that("psislw returns expected results", {
  psis <- psislw(x[, 1])
  lw <- psis$lw_smooth
  expect_equal(length(psis), 2L)
  expect_equal(nrow(lw), nrow(x))
  expect_equal(lw[1], -5.6655489517740527106)
  expect_equal(lw[50], -5.188442371693668953)
  expect_equal(range(lw), c(-7.4142421808626526314, -2.6902215137943321643))
  expect_equal(psis$pareto_k, 0.17364505906017813075)
})

test_that("psislw handles special cases, throws appropriate errors/warnings", {
  expect_warning(psis <- psislw(x[, 1], wcp = 0.01),
                 regexp = "Too few tail samples")
  expect_true(is.infinite(psis$pareto_k))
  expect_error(psislw(wcp = 0.2),
               regexp = "'lw' or 'llfun' and 'llargs' must be specified")
})
