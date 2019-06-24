library(loo)

context("psis_approximate_posterior")

load(test_path("test_data_psis_approximate_posterior.rda"))


test_that("paloo.array works as paloo.matrix", {
  skip_if_not_installed("checkmate")
  log_p <- test_data_psis_approximate_posterior$laplace_independent$log_p
  log_g <- test_data_psis_approximate_posterior$laplace_independent$log_q
  ll <- test_data_psis_approximate_posterior$laplace_independent$log_liks

  # Create array with two "chains"
  log_p_mat <- matrix(log_p, nrow = 500, ncol = 2)
  log_g_mat <- matrix(log_g, nrow = 500, ncol = 2)
  ll_array <- array(0, dim = c(500, 2 ,ncol(ll)))
  ll_array[,1,] <- ll[1:500,]
  ll_array[,2,] <- ll[501:1000,]

  # Assert that they are ok
  expect_equivalent(ll_array[1:2,1,1:2], ll[1:2,1:2])
  expect_equivalent(ll_array[1:2,2,1:2], ll[501:502,1:2])

  # Compute aploo
  expect_silent(aploo1 <- aploo.matrix(x = ll, log_p = log_p, log_g = log_g))
  expect_silent(aploo2 <- aploo.array(x = ll_array, log_p = log_p_mat, log_g = log_g_mat))
  expect_silent(aploo1b <- loo.matrix(x = ll, r_eff = rep(1, 10)))

  # Check equivalence
  expect_equal(aploo1$estimates, aploo2$estimates)
  expect_failure(expect_equal(aploo1b$estimates, aploo2$estimates))
})


