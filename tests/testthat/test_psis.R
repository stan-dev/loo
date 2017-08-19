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
  psis_with_fn <- psis(llfun, args = nlist(data, draws, N, S))

  # Check that we get same answer if using log-likelihood matrix
  log_lik_mat <- sapply(1:N, function(i) llfun(i, data[i,, drop=FALSE], draws))
  psis_with_mat <- loo(log_lik_mat)
  expect_equal(psis_with_mat, psis_with_fn)
})

