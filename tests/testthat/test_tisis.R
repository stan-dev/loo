library(loo)
options(mc.cores=1)
options(loo.cores=NULL)
set.seed(123)

context("tis and is")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = dim(LLarr)[1])
r_eff_arr <- relative_eff(exp(LLarr))
r_eff_vec <- relative_eff(exp(LLvec), chain_id = chain_id)
psis1 <- psis(log_ratios = -LLarr, r_eff = r_eff_arr)
tis1 <- loo:::importance_sampling(log_ratios = -LLarr, r_eff = r_eff_arr, is_method = "TIS")
is1 <- loo:::importance_sampling(log_ratios = -LLarr, r_eff = r_eff_arr, is_method = "IS")

test_that("tis and is runs", {
  LLvec[1] <- -10
  expect_silent(tis1 <- loo:::importance_sampling(log_ratios = -LLvec, r_eff = r_eff_vec, is_method = "TIS"))
  expect_silent(is1 <- loo:::importance_sampling(log_ratios = -LLvec, r_eff = r_eff_vec, is_method = "IS"))
  expect_failure(expect_equal(tis1$log_weights, is1$log_weights))
  expect_failure(expect_equal(tis1$log_weights, psis1$log_weights))
})

test_that("loo:::importance_sampling returns object with correct structure for TIS/IS", {
  skip("Should be fixed")
  expect_true(is.psis(tis1))
  expect_true(is.psis(is1))

  expect_false(is.loo(tis1))
  expect_false(is.loo(is1))

  expect_false(is.psis_loo(tis1))
  expect_false(is.psis_loo(is1))

  expect_named(tis1, c("log_weights", "diagnostics"))
  expect_named(is1, c("log_weights", "diagnostics"))

  expect_named(tis1$diagnostics, c("pareto_k", "n_eff"))
  expect_named(is1$diagnostics, c("pareto_k", "n_eff"))

  expect_equal(dim(tis1), dim(LLmat))
  expect_equal(dim(is1), dim(LLmat))

  expect_length(tis1$diagnostics$pareto_k, dim(psis1)[2])
  expect_length(is1$diagnostics$pareto_k, dim(psis1)[2])

  expect_length(tis1$diagnostics$n_eff, dim(psis1)[2])
  expect_length(is1$diagnostics$n_eff, dim(psis1)[2])

  expect_equal(attr(psis1, "is_method")[1], "PSIS")
  expect_equal(attr(tis1, "is_method")[1], "TIS")
  expect_equal(attr(is1, "is_method")[1], "IS")
})


test_that("psis methods give same results", {
  tis2 <- suppressWarnings(loo:::importance_sampling(-LLmat, r_eff = r_eff_arr, is_method = "TIS"))
  expect_identical(tis1, tis2)

  tisvec <- suppressWarnings(loo:::importance_sampling(-LLvec, r_eff = r_eff_vec, is_method = "TIS"))
  tismat <- suppressWarnings(loo:::importance_sampling(-LLmat[, 1], r_eff = r_eff_vec, is_method = "TIS"))
  expect_identical(tisvec, tismat)

  is2 <- suppressWarnings(loo:::importance_sampling(-LLmat, r_eff = r_eff_arr, is_method = "IS"))
  expect_identical(is1, is2)

  isvec <- suppressWarnings(loo:::importance_sampling(-LLvec, r_eff = r_eff_vec, is_method = "IS"))
  ismat <- suppressWarnings(loo:::importance_sampling(-LLmat[, 1], r_eff = r_eff_vec, is_method = "IS"))
  expect_identical(isvec, ismat)
})



test_that("psis throws correct errors and warnings", {
  # r_eff=NULL warnings
  expect_warning(loo:::importance_sampling(-LLarr, is_method = "TIS"), "Relative effective sample sizes")
  expect_warning(loo:::importance_sampling(-LLmat, is_method = "TIS"), "Relative effective sample sizes")
  expect_warning(loo:::importance_sampling(-LLmat[, 1], is_method = "TIS"), "Relative effective sample sizes")

  # r_eff=NA disables warnings
  expect_silent(loo:::importance_sampling(-LLarr, r_eff = NA, is_method = "TIS"))
  expect_silent(loo:::importance_sampling(-LLmat, r_eff = NA, is_method = "TIS"))
  expect_silent(loo:::importance_sampling(-LLmat[,1], r_eff = NA, is_method = "TIS"))

  # r_eff=NULL and r_eff=NA give same answer
  expect_equal(
    suppressWarnings(loo:::importance_sampling(-LLarr, is_method = "TIS")),
    loo:::importance_sampling(-LLarr, r_eff = NA, is_method = "TIS")
  )

  # r_eff wrong length is error
  expect_error(loo:::importance_sampling(-LLarr, r_eff = r_eff_arr[-1], is_method = "TIS"), "one value per observation")

  # r_eff has some NA values causes error
  r_eff_arr[2] <- NA
  expect_error(loo:::importance_sampling(-LLarr, r_eff = r_eff_arr, is_method = "TIS"), "mix NA and not NA values")

  # no NAs or non-finite values allowed
  LLmat[1,1] <- NA
  expect_error(loo:::importance_sampling(-LLmat, is_method = "TIS"), "NAs not allowed in input")

  LLmat[1,1] <- 1
  LLmat[10, 2] <- Inf
  expect_error(loo:::importance_sampling(-LLmat, is_method = "TIS"), "All input values must be finite")

  # no lists allowed
  expect_error(expect_warning(loo:::importance_sampling(as.list(-LLvec), is_method = "TIS"), "List not allowed as input"))

  # if array, must be 3-D array
  dim(LLarr) <- c(2, 250, 2, 32)
  expect_error(
    loo:::importance_sampling(-LLarr, is_method = "TIS"),
    "length(dim(log_ratios)) == 3 is not TRUE",
    fixed = TRUE
  )
})


test_that("explict test of values for IS and TIS", {
  lw <- 1:16
  # With TIS values greater than 0.5*log(length(lw)) i.e. 1.386294, should be truncated to 10
  expect_silent(tis_true <- loo:::importance_sampling(log_ratios = lw, r_eff = NA, is_method = "TIS"))
  expect_equal(as.vector(weights(tis_true, log = TRUE, normalize = FALSE)),
               c(-14.0723, -13.0723, -12.0723, -11.0723, -10.0723, -9.0723, -8.0723, -7.0723, -6.0723, -5.0723, -4.0723, -3.0723, -2.0723, -1.0723, -0.0723, 0.), tol = 0.001)
  expect_silent(is_true <- loo:::importance_sampling(log_ratios = lw, r_eff = NA, is_method = "IS"))
  expect_equal(as.vector(weights(is_true, log = TRUE, normalize = FALSE)),
               c(-15:0), tol = 0.00001)

  lw <- c(0.7609420, 1.3894140, 0.4158346, 2.5307927, 4.3379119, 2.4159240, 2.2462172, 0.8057697, 0.9333107, 1.5599302)

  expect_silent(tis_true <- loo:::importance_sampling(log_ratios = lw, r_eff = NA, is_method = "TIS"))
  expect_equal(as.vector(weights(tis_true, log = TRUE, normalize = FALSE)),
               c(-2.931, -2.303, -3.276, -1.161,  0, -1.276, -1.446, -2.886, -2.759, -2.132),
               tol = 0.001)
  expect_silent(is_true <- loo:::importance_sampling(log_ratios = lw, r_eff = NA, is_method = "IS"))
  expect_equal(as.vector(weights(is_true, log = TRUE, normalize = FALSE)),
               lw-max(lw), tol = 0.00001)
})


