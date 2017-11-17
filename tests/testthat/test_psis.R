library(loo)
options(loo.cores = 2)
set.seed(123)

context("psis")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = dim(LLarr)[1])
r_eff_arr <- relative_eff(exp(LLarr))
r_eff_vec <- relative_eff(exp(LLvec), chain_id = chain_id)
psis1 <- psis(log_ratios = -LLarr, r_eff = r_eff_arr)


test_that("psis returns object with correct structure", {
  expect_true(is.psis(psis1))
  expect_false(is.loo(psis1))
  expect_false(is.psis_loo(psis1))

  expect_named(psis1, c("log_weights", "diagnostics"))
  expect_named(psis1$diagnostics, c("pareto_k", "n_eff"))
  expect_equal(dim(psis1), dim(LLmat))
  expect_length(psis1$diagnostics$pareto_k, dim(psis1)[2])
  expect_length(psis1$diagnostics$n_eff, dim(psis1)[2])
})


test_that("psis methods give same results", {
  psis2 <- suppressWarnings(psis(-LLmat, r_eff = r_eff_arr))
  expect_identical(psis1, psis2)

  psisvec <- suppressWarnings(psis(-LLvec, r_eff = r_eff_vec))
  psismat <- suppressWarnings(psis(-LLmat[, 1], r_eff = r_eff_vec))
  expect_identical(psisvec, psismat)
})

test_that("psis throws correct errors", {
  # no NAs or non-finite values allowed
  LLmat[1,1] <- NA
  expect_error(psis(-LLmat), "NAs not allowed in input")

  LLmat[1,1] <- 1
  LLmat[10, 2] <- Inf
  expect_error(psis(-LLmat), "All input values must be finite")

  # no lists allowed
  expect_error(psis(as.list(-LLvec)), "List not allowed as input")

  # if array, must be 3-D array
  dim(LLarr) <- c(2, 250, 2, 32)
  expect_error(
    psis(-LLarr),
    "length(dim(log_ratios)) == 3 is not TRUE",
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


test_that("psis_n_eff methods works properly", {
  w <- weights(psis1, normalize = TRUE, log = FALSE)
  expect_equal(psis_n_eff.default(w[, 1], r_eff = 1), 1 / sum(w[, 1]^2))
  expect_equal(psis_n_eff.default(w[, 1], r_eff = 2), 2 / sum(w[, 1]^2))
  expect_equal(
    psis_n_eff.default(w[, 1], r_eff = 2),
    psis_n_eff.matrix(w, r_eff = rep(2, ncol(w)))[1]
  )
  expect_warning(psis_n_eff.default(w[, 1]), "not adjusted based on MCMC n_eff")
  expect_warning(psis_n_eff.matrix(w), "not adjusted based on MCMC n_eff")
})

test_that("default relative_eff method works properly", {
  expect_equal(relative_eff.default(LLmat[, 1], chain_id),
               mcmc_n_eff(LLarr[, , 1]) / 1000)
})

