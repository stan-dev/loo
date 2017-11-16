library(loo)
options(loo.cores = 1)
set.seed(123)

context("loo and waic")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = nrow(LLarr))

r_eff_arr <- relative_eff(exp(LLarr))
r_eff_mat <- relative_eff(exp(LLmat), chain_id = chain_id)

l1 <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr))
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
  expect_equal(dim(w1), dim(LLmat))
})

test_that("loo returns object with correct structure", {
  expect_false(is.waic(l1))
  expect_true(is.loo(l1))
  expect_true(is.psis_loo(l1))
  expect_named(l1, c("estimates", "pointwise", "diagnostics", "psis_object"))
  expect_named(l1$diagnostics, c("pareto_k", "n_eff"))
  expect_equal(dimnames(l1$estimates)[[1]], c("elpd_loo", "p_loo", "looic"))
  expect_equal(dimnames(l1$estimates)[[2]], c("Estimate", "SE"))
  expect_equal(colnames(l1$pointwise), c("elpd_loo", "mcse_elpd_loo", "p_loo", "looic"))
  expect_equal(dim(l1), dim(LLmat))
})

test_that("loo.array and loo.matrix give same result", {
  l2 <- suppressWarnings(loo(LLmat, r_eff = r_eff_mat))
  expect_identical(l1$estimates, l2$estimates)
  expect_identical(l1$diagnostics, l2$diagnostics)

  # the mcse_elpd_loo columns won't be identical because we use sampling
  expect_identical(l1$pointwise[, -2], l2$pointwise[, -2])
  expect_equal(l1$pointwise[, 2], l2$pointwise[, 2], tol = 0.005)
})

test_that("waic.array and waic.matrix give same result", {
  w2 <- suppressWarnings(waic(LLmat))
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
  draws <- as.matrix(rbeta(S, a, b))
  data <- data.frame(y,K)
  llfun <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }
  loo_with_fn <- loo(llfun, data = data, draws = draws, r_eff = rep(1, nrow(data)))
  waic_with_fn <- waic(llfun, data = data, draws = draws)

  # Check that we get same answer if using log-likelihood matrix
  log_lik_mat <- sapply(1:N, function(i) llfun(data[i,, drop=FALSE], draws))
  loo_with_mat <- loo(log_lik_mat, r_eff = rep(1, ncol(log_lik_mat)))
  waic_with_mat <- waic(log_lik_mat)

  expect_equal(waic_with_mat, waic_with_fn)
  expect_identical(loo_with_mat$estimates, loo_with_fn$estimates)
  expect_identical(loo_with_mat$diagnostics, loo_with_fn$diagnostics)
})

