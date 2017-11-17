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

loo1 <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr))
waic1 <- suppressWarnings(waic(LLarr))


test_that("loo and waic results haven't changed", {
  expect_equal_to_reference(loo1, "loo.rds")
  expect_equal_to_reference(waic1, "waic.rds")
})

test_that("loo with cores=1 and cores=2 gives same results", {
  loo2 <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr, cores = 2))
  expect_equal(loo1$estimates, loo2$estimates)
})

test_that("waic returns object with correct structure", {
  expect_true(is.waic(waic1))
  expect_true(is.loo(waic1))
  expect_false(is.psis_loo(waic1))
  expect_named(waic1, c("estimates", "pointwise"))
  est_names <- dimnames(waic1$estimates)
  expect_equal(est_names[[1]], c("elpd_waic", "p_waic", "waic"))
  expect_equal(est_names[[2]], c("Estimate", "SE"))
  expect_equal(colnames(waic1$pointwise), est_names[[1]])
  expect_equal(dim(waic1), dim(LLmat))
})

test_that("loo returns object with correct structure", {
  expect_false(is.waic(loo1))
  expect_true(is.loo(loo1))
  expect_true(is.psis_loo(loo1))
  expect_named(loo1, c("estimates", "pointwise", "diagnostics", "psis_object"))
  expect_named(loo1$diagnostics, c("pareto_k", "n_eff"))
  expect_equal(dimnames(loo1$estimates)[[1]], c("elpd_loo", "p_loo", "looic"))
  expect_equal(dimnames(loo1$estimates)[[2]], c("Estimate", "SE"))
  expect_equal(colnames(loo1$pointwise), c("elpd_loo", "mcse_elpd_loo", "p_loo", "looic"))
  expect_equal(dim(loo1), dim(LLmat))
})

test_that("loo.array and loo.matrix give same result", {
  l2 <- suppressWarnings(loo(LLmat, r_eff = r_eff_mat))
  expect_identical(loo1$estimates, l2$estimates)
  expect_identical(loo1$diagnostics, l2$diagnostics)

  # the mcse_elpd_loo columns won't be identical because we use sampling
  expect_identical(loo1$pointwise[, -2], l2$pointwise[, -2])
  expect_equal(loo1$pointwise[, 2], l2$pointwise[, 2], tol = 0.005)
})

test_that("waic.array and waic.matrix give same result", {
  waic2 <- suppressWarnings(waic(LLmat))
  expect_identical(waic1, waic2)
})

test_that("loo and waic error with vector input", {
  expect_error(loo(LLvec), regexp = "no applicable method")
  expect_error(waic(LLvec), regexp = "no applicable method")
})



# fake data and posterior draws for testing function methods
N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
p <- rbeta(1, a0, b0)
y <- rbinom(N, size = K, prob = p)
a <- a0 + sum(y); b <- b0 + N * K - sum(y)
draws <- as.matrix(rbeta(S, a, b))
data <- data.frame(y,K)
llfun <- function(data_i, draws) {
  dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
}
log_lik_mat <- sapply(1:N, function(i) llfun(data[i,, drop=FALSE], draws))

waic_with_fn <- waic(llfun, data = data, draws = draws)
waic_with_mat <- waic(log_lik_mat)

loo_with_fn <- loo(llfun, data = data, draws = draws,
                   r_eff = rep(1, nrow(data)))
loo_with_mat <- loo(log_lik_mat, r_eff = rep(1, ncol(log_lik_mat)),
                    save_psis = TRUE)

test_that("loo_i results match loo results for ith data point", {
  loo_i_val <- loo_i(i = 2, llfun = llfun, data = data, draws = draws)
  expect_equal(loo_i_val$pointwise[, "elpd_loo"], loo_with_fn$pointwise[2, "elpd_loo"])
  expect_equal(loo_i_val$pointwise[, "p_loo"], loo_with_fn$pointwise[2, "p_loo"])
  expect_equal(loo_i_val$diagnostics$pareto_k, loo_with_fn$diagnostics$pareto_k[2])
  expect_equal(loo_i_val$diagnostics$n_eff, loo_with_fn$diagnostics$n_eff[2])
})

test_that("function and matrix methods return same result", {
  expect_equal(waic_with_mat, waic_with_fn)
  expect_identical(loo_with_mat$estimates, loo_with_fn$estimates)
  expect_identical(loo_with_mat$diagnostics, loo_with_fn$diagnostics)
  expect_identical(dim(loo_with_mat), dim(loo_with_fn))
})

test_that("save_psis option to loo.function makes correct psis object", {
  loo_with_fn2 <- loo(llfun, data = data, draws = draws,
                      r_eff = rep(1, nrow(data)), save_psis = TRUE)
  expect_identical(loo_with_fn2$psis_object, loo_with_mat$psis_object)
})

