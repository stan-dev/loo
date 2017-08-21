library(loo)
options(loo.cores = 2)
set.seed(123)

context("psis")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = dim(LLarr)[1])

psis1 <- suppressWarnings(psis(LLarr))

test_that("psis returns object with correct structure", {
  expect_true(is.psis(psis1))
  expect_false(is.loo(psis1))
  expect_false(is.psis_loo(psis1))

  expect_named(psis1, c("log_weights", "pareto_k", "n_eff"))
  expect_equal(dim(psis1$log_weights), dim(LLmat))
  expect_equal(length(psis1$pareto_k), ncol(psis1$log_weights))
  expect_equal(length(psis1$n_eff), ncol(psis1$log_weights))
})

test_that("psis methods give same results", {
  psis2 <- suppressWarnings(psis(LLmat, chain_id))
  expect_identical(psis1, psis2)

  psisvec <- suppressWarnings(psis(LLvec, chain_id))
  psismat <- suppressWarnings(psis(LLmat[, 1], chain_id))
  expect_identical(psisvec, psismat)
})

test_that("psis throws correct errors", {
  # vector and matrix methods need chain_id
  expect_error(psis(LLvec), 'argument "chain_id" is missing')
  expect_error(psis(LLmat), 'argument "chain_id" is missing')

  # no NAs allowed
  LLmat[1,1] <- NA
  expect_error(
    psis(LLmat, chain_id),
    "!anyNA(x) is not TRUE",
    fixed = TRUE
  )

  # if array, must be 3-D array
  dim(LLarr) <- c(2, 250, 2, 32)
  expect_error(
    psis(LLarr),
    "length(dim(x)) == 3 is not TRUE",
    fixed = TRUE
  )
})


test_that("weights method returns correct output", {
  # default arguments
  expect_identical(weights(psis1), weights(psis1, normalize = TRUE, log = TRUE))

  # unnormalized log-weights same as in psis object
  expect_equal(psis1$log_weights, weights(psis1, normalize = FALSE))

  # normalized weights sum to 1
  expect_equal(
    colSums(weights(psis1, normalize = TRUE, log = FALSE)),
    rep(1, ncol(psis1$log_weights))
  )
})


test_that("default psis_n_eff method works properly", {
  w <- weights(psis1, normalize = TRUE, log = FALSE)
  expect_equal(psis_n_eff.default(w[, 1], rel_n_eff = 1), 1 / sum(w[, 1]^2))
  expect_equal(psis_n_eff.default(w[, 1], rel_n_eff = 2), 2 / sum(w[, 1]^2))
  expect_warning(psis_n_eff.default(w[, 1]), "not adjusted based on MCMC n_eff")
})

test_that("default relative_n_eff method works properly", {
  expect_equal(relative_n_eff.default(LLmat[, 1], chain_id),
               mcmc_n_eff(LLarr[, , 1]) / 1000)
})

