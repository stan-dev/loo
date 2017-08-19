library(loo)
options(loo.cores = 1)
set.seed(123)

context("loo and waic")

LLarr <- source(test_path("LL_array_data.R"))$value
LLmat <- llarray_to_matrix(LLarr)
LLvec <- LLmat[, 1]

test_that("loo.array and loo.matrix give same result", {
  l1 <- suppressWarnings(loo(LLarr))
  l2 <- suppressWarnings(loo(LLmat, chain_id = rep(1:2, each = 50)))
  expect_identical(l1, l2)
})

test_that("waic.array and loo.matrix give same result", {
  w1 <- suppressWarnings(waic(LLarr))
  w2 <- suppressWarnings(waic(LLmat, chain_id = rep(1:2, each = 50)))
  expect_identical(w1, w2)
})


test_that("loo and waic error with vector input", {
  expect_error(loo(LLvec), regexp = "no applicable method")
  expect_error(waic(LLvec), regexp = "no applicable method")
})

test_that("function and matrix methods return same result", {
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
  # loo_with_fn <- loo(llfun, args = nlist(data, draws, N, S))
  waic_with_fn <- waic(llfun, args = nlist(data, draws, N, S))

  # Check that we get same answer if using log-likelihood matrix
  log_lik_mat <- sapply(1:N, function(i) llfun(i, data[i,, drop=FALSE], draws))
  # loo_with_mat <- loo(log_lik_mat)
  waic_with_mat <- waic(log_lik_mat)
  # expect_equal(loo_with_mat, loo_with_fn)
  expect_equal(waic_with_mat, waic_with_fn)
})

