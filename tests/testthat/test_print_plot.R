library(loo)
options(loo.cores = 1)

set.seed(123)
x1 <- waic(matrix(rnorm(5000), 100, 50))
x2 <- loo(matrix(rnorm(5000), 100, 50))

# test loo and waic -------------------------------------------------------
context("print and plot")
test_that("plot.loo throws appropriate errors", {
  expect_error(plot(x1), regexp = "No Pareto k values")
})
test_that("print.loo issues appropriate warnings for waic",{
  llmsg <- "Computed from 100 by 50 log-likelihood matrix"
  if (any(x1$pointwise[, "p_waic"] > 0.4)) expect_warning(print(x1))
  else expect_output(print(x1), regexp = llmsg)
  x1$pointwise[, "p_waic"] <- runif(50, 0, .3)
  expect_output(print(x1), regexp = llmsg)
  x1$pointwise[1:5, "p_waic"] <- 0.5
  expect_warning(print(x1), regexp = "p_waic estimates greater than 0.4")
})
test_that("print.loo issues appropriate warnings for loo",{
  if (any(x2$pareto_k > 0.5)) expect_warning(print(x2))
  else expect_output(print(x2), regexp = "All Pareto k estimates OK")
  x2$pareto_k <- runif(50, 0, .5)
  expect_output(print(x2), regexp = "All Pareto k estimates OK")
  x2$pareto_k[1] <- 0.51
  expect_warning(print(x2), regexp = "between 0.5 and 1")
  x2$pareto_k[1] <- 1.1
  expect_warning(print(x2), regexp = "greater than 1")
})
