library(loo)
options(loo.cores = 1)

set.seed(123)
x1 <- waic(matrix(rnorm(5000), 100, 50))
x2 <- loo(matrix(rnorm(5000), 100, 50))

context("print and plot")
test_that("plot.loo throws appropriate errors", {
  expect_error(plot(x1), regexp = "No Pareto k estimates found")
})
test_that("plot.loo doesn't error", {
  # this doesn't actually check that plots are the same, but just serves
  # to check that they don't throw errors
  expect_equal(plot(x2), loo:::plot_k(x2$pareto_k))
  expect_equal(plot(x2), loo:::plot_k(x2$pareto_k))
  expect_equal(plot(x2, label_points = TRUE),
               loo:::plot_k(x2$pareto_k, label_points = TRUE))
  x3 <- x2
  x3$pareto_k[1:5] <- Inf
  expect_warning(plot(x3), regexp = "estimates are Inf/NA/NaN and not plotted.")
})
test_that("print.loo issues appropriate warnings for waic",{
  llmsg <- "Computed from 100 by 50 log-likelihood matrix"
  if (any(x1$pointwise[, "p_waic"] > 0.4))
    expect_warning(print(x1), "p_waic estimates greater than 0.4")

  expect_output(print(x1), regexp = llmsg)

  x1$pointwise[, "p_waic"] <- runif(50, 0, .3)
  expect_output(print(x1), regexp = llmsg)
  x1$pointwise[1:5, "p_waic"] <- 0.5
  expect_warning(print(x1), regexp = "p_waic estimates greater than 0.4")
})
test_that("print.loo works correctly for loo",{
  if (any(x2$pareto_k > 0.7)) {
    expect_output(print(x2), "Pareto k diagnostic values")
  } else if (any(x2$pareto_k > 0.5)) {
    expect_output(print(x2), regexp = "All Pareto k estimates are ok")
  } else {
    expect_output(print(x2), regexp = "All Pareto k estimates are good")
  }

  x2$pareto_k <- runif(50, 0, .49)
  expect_output(print(x2), regexp = "All Pareto k estimates are good")
  x2$pareto_k[1] <- 0.71
  expect_output(print(x2), regexp = "Pareto k diagnostic values")
  x2$pareto_k[1] <- 1.1
  expect_output(print(x2), regexp = "Pareto k diagnostic values")
})
