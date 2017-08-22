library(loo)
options(loo.cores = 1)
set.seed(123)

context("loo and waic")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = nrow(LLarr))

l1 <- suppressWarnings(loo(LLarr))
w1 <- suppressWarnings(waic(LLarr))

test_that("loo and waic results haven't changed", {
  expect_equal_to_reference(l1, "loo.rds")
  expect_equal_to_reference(w1, "waic.rds")
})

test_that("waic returns object with correct structure", {
  expect_true(is.waic(w1))
  expect_true(is.loo(w1))
  expect_false(is.psis_loo(w1))
  expect_named(w1, c("estimates", "pointwise"))
  est_names <- dimnames(w1$estimates)
  expect_equal(est_names[[1]], c("elpd_waic", "p_waic", "waic"))
  expect_equal(est_names[[2]], c("Estimate", "SE"))
  expect_equal(colnames(w1$pointwise), est_names[[1]])
})

test_that("loo returns object with correct structure", {
  expect_false(is.waic(l1))
  expect_true(is.loo(l1))
  expect_true(is.psis_loo(l1))
  expect_named(l1, c("estimates", "pointwise", "diagnostics"))
  expect_named(l1$diagnostics, c("pareto_k", "n_eff"))
  est_names <- dimnames(l1$estimates)
  expect_equal(est_names[[1]], c("elpd_loo", "p_loo", "looic"))
  expect_equal(est_names[[2]], c("Estimate", "SE"))
  expect_equal(colnames(l1$pointwise), est_names[[1]])
})

test_that("loo.array and loo.matrix give same result", {
  l2 <- suppressWarnings(loo(LLmat, chain_id = chain_id))
  expect_identical(l1, l2)
})

test_that("waic.array and waic.matrix give same result", {
  w2 <- suppressWarnings(waic(LLmat, chain_id = chain_id))
  expect_identical(w1, w2)
})

test_that("loo and waic error with vector input", {
  expect_error(loo(LLvec), regexp = "no applicable method")
  expect_error(waic(LLvec), regexp = "no applicable method")
})

test_that("waic function and matrix methods return same result", {
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

