library(loo)
SW <- suppressWarnings

context("psislw")

set.seed(123)
x <- matrix(rnorm(5000), 100, 50)

expect_deprecated <- function(object) {
  testthat::expect_warning(object, "deprecated", ignore.case = TRUE)
}

test_that("psislw throws deprecation warning", {
  expect_deprecated(psislw(x[, 1]))
})


test_that("psislw handles special cases, throws appropriate errors/warnings", {
  expect_warning(
    psis <- psislw(x[, 1], wcp = 0.01),
    regexp = "All tail values are the same. Weights are truncated but not smoothed"
  )
  expect_true(is.infinite(psis$pareto_k))

  expect_warning(
    psislw(x[, 1], wcp = 0.01),
    regexp = "Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details",
    fixed = TRUE
  )

  expect_error(
    expect_deprecated(psislw(wcp = 0.2)),
    regexp = "'lw' or 'llfun' and 'llargs' must be specified"
  )
})

test_that("psislw returns expected results", {
  psis <- SW(psislw(x[, 1]))
  lw <- psis$lw_smooth
  expect_equal(length(psis), 2L)
  expect_equal(nrow(lw), nrow(x))
  expect_equal(lw[1], -5.6655489517740527106)
  expect_equal(lw[50], -5.188442371693668953)
  expect_equal(range(lw), c(-7.4142421808626526314, -2.6902215137943321643))
  expect_equal(psis$pareto_k, 0.17364505906017813075)
})

test_that("psislw function and matrix methods return same result", {
  set.seed(024)

  # fake data and posterior draws
  N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
  p <- rbeta(1, a0, b0)
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  draws <- rbeta(S, a, b)
  data <- data.frame(y,K)
  llfun <- function(i, data, draws) {
    dbinom(data$y, size = data$K, prob = draws, log = TRUE)
  }
  psislw_with_fn <- SW(psislw(llfun = llfun, llargs = nlist(data, draws, N, S)))

  # Check that we get same answer if using log-likelihood matrix
  ll <- sapply(1:N, function(i) llfun(i, data[i,, drop=FALSE], draws))
  psislw_with_mat <- SW(psislw(-ll))
  expect_equal(psislw_with_fn, psislw_with_mat)
})

test_that("psislw_warnings helper works properly", {
  k <- c(0, 0.1, 0.55, 0.75)
  expect_silent(psislw_warnings(k[1:2]))
  expect_warning(psislw_warnings(k[1:3]),
                 "Some Pareto k diagnostic values are slightly high")
  expect_warning(psislw_warnings(k),
                 "Some Pareto k diagnostic values are too high")
})
