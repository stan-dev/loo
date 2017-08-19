library(loo)
context("helper functions")

LLarr <- source(test_path("LL_array_data.R"))$value

test_that("reshaping functions result in correct dimensions", {
  LLmat <- llarray_to_matrix(LLarr)
  expect_equal(dim(LLmat), c(100, 32))

  LLarr2 <- llmatrix_to_array(LLmat, chain_id = rep(1:2, each = 50))
  expect_identical(LLarr, LLarr2)
})

test_that("logColMeansExp(x) = log(colMeans(exp(x))) ", {
  x <- matrix(rnorm(100), 20, 5)
  expect_equal(logColMeansExp(x), log(colMeans(exp(x))))
})

test_that("nlist works", {
  a <- 1; b <- 2; c <- 3;
  nlist_val <- list(nlist(a, b, c), nlist(a, b, c = "tornado"))
  nlist_ans <- list(list(a = 1, b = 2, c = 3), list(a = 1, b = 2, c = "tornado"))
  expect_equal(nlist_val, nlist_ans)
})
