library(loo)
options(mc.cores = 1)
set.seed(123)

context("relative_eff methods")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()

test_that("relative_eff results haven't changed", {
  expect_equal_to_reference(relative_eff(exp(LLarr)), "reference-results/relative_eff.rds")
})

test_that("relative_eff is equal to ESS / S", {
  dims <- dim(LLarr)
  ess <- r_eff <- rep(NA, dims[3])
  for (j in 1:dims[3]) {
    r_eff[j] <- relative_eff(LLarr[,,1, drop=FALSE])
    ess[j] <- ess_rfun(LLarr[,,1])
  }
  S <- prod(dim(LLarr)[1:2])
  expect_equal(r_eff, ess / S)
})

test_that("relative_eff array and matrix methods return identical output", {
  r_eff_arr <- relative_eff(exp(LLarr))
  r_eff_mat <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = nrow(LLarr)))
  expect_identical(r_eff_arr, r_eff_mat)
})

test_that("relative_eff matrix and function methods return identical output", {
  source(test_path("data-for-tests/function_method_stuff.R"))
  chain <- rep(1, nrow(draws))
  r_eff_mat <- relative_eff(llmat_from_fn, chain_id = chain)
  r_eff_fn <- relative_eff(llfun, chain_id = chain, data = data, draws = draws, cores = 1)
  expect_identical(r_eff_mat, r_eff_fn)
})

test_that("relative_eff with multiple cores runs", {
  skip_on_cran()
  source(test_path("data-for-tests/function_method_stuff.R"))
  dim(llmat_from_fn) <- c(nrow(llmat_from_fn), 1, ncol(llmat_from_fn))
  r_eff_arr <- relative_eff(llmat_from_fn, cores = 2)
  r_eff_fn <-
    relative_eff(
      llfun,
      chain_id = rep(1, nrow(draws)),
      data = data,
      draws = draws,
      cores = 2
    )
  expect_identical(r_eff_arr, r_eff_fn)
})
