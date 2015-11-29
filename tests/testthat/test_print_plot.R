library(loo)
options(loo.cores = 1)

# test loo and waic -------------------------------------------------------
context("print and plot")
test_that("plot.loo throws appropriate errors", {
  set.seed(123)
  x1 <- waic(matrix(rnorm(5000), 100, 50))
  x2 <- loo(matrix(rnorm(5000), 100, 50))
  expect_error(plot(x1), regexp = "No Pareto k values")

  if (any(x2$pareto_k > 0.5)) expect_warning(print(x2))
  else expect_output(print(x2), regexp = "All Pareto k estimates OK")

  x2$pareto_k <- runif(50, 0, .5)
  expect_output(print(x2), regexp = "All Pareto k estimates OK")
  x2$pareto_k[1] <- 0.51
  expect_warning(print(x2), regexp = "between 0.5 and 1")
  x2$pareto_k[1] <- 1.1
  expect_warning(print(x2), regexp = "greater than 1")
})
