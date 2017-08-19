library(loo)
options(loo.cores = 2)
set.seed(123)

context("psis")

LLarr <- source(test_path("LL_array_data.R"))$value
LLmat <- llarray_to_matrix(LLarr)
LLvec <- LLmat[, 1]

psis1 <- suppressWarnings(psis(LLarr))

test_that("psis methods give same results", {
  psis2 <- suppressWarnings(psis(LLmat, chain_id = rep(1:2, each = 50)))
  expect_identical(psis1, psis2)

  psisvec <- suppressWarnings(psis(LLvec, chain_id = rep(1:2, each = 50)))
  psismat <- suppressWarnings(psis(LLmat[, 1], chain_id = rep(1:2, each = 50)))
  expect_identical(psisvec, psismat)
})

test_that("psis throws correct errors", {
  # vector and matrix methods need chain_id
  expect_error(psis(LLvec), 'argument "chain_id" is missing')
  expect_error(psis(LLmat), 'argument "chain_id" is missing')

  # no NAs allowed
  LLmat[1,1] <- NA
  expect_error(
    psis(LLmat, chain_id = rep(1:2, each = 50)),
    "!anyNA(x) is not TRUE",
    fixed = TRUE
  )

  # if array, must be 3-D array
  dim(LLarr) <- c(2, 25, 2, 32)
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

