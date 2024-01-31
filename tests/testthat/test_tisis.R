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
tis1 <- tis(log_ratios = -LLarr, r_eff = r_eff_arr)
is1 <- sis(log_ratios = -LLarr, r_eff = r_eff_arr)


test_that("tis and is runs", {
  LLvec[1] <- -10
  expect_silent(tis1 <- tis(log_ratios = -LLvec, r_eff = r_eff_vec))
  expect_silent(is1 <- sis(log_ratios = -LLvec, r_eff = r_eff_vec))
  expect_failure(expect_equal(tis1$log_weights, is1$log_weights))
  expect_failure(expect_equal(tis1$log_weights, psis1$log_weights))
})

test_that("tis() and sis() returns object with correct structure for tis/sis", {

  expect_false(is.psis(tis1))
  expect_false(is.psis(is1))
  expect_true(is.tis(tis1))
  expect_false(is.tis(is1))
  expect_false(is.sis(tis1))
  expect_true(is.sis(is1))

  expect_false(is.loo(tis1))
  expect_false(is.loo(is1))

  expect_false(is.psis_loo(tis1))
  expect_false(is.psis_loo(is1))

  expect_named(tis1, c("log_weights", "diagnostics"))
  expect_named(is1, c("log_weights", "diagnostics"))

  expect_named(tis1$diagnostics, c("pareto_k", "n_eff", "r_eff"))
  expect_named(is1$diagnostics, c("pareto_k", "n_eff", "r_eff"))

  expect_equal(dim(tis1), dim(LLmat))
  expect_equal(dim(is1), dim(LLmat))

  expect_length(tis1$diagnostics$pareto_k, dim(psis1)[2])
  expect_length(is1$diagnostics$pareto_k, dim(psis1)[2])

  expect_length(tis1$diagnostics$n_eff, dim(psis1)[2])
  expect_length(is1$diagnostics$n_eff, dim(psis1)[2])

  expect_equal(attr(psis1, "method")[1], "psis")
  expect_equal(attr(tis1, "method")[1], "tis")
  expect_equal(attr(is1, "method")[1], "sis")
})


test_that("psis methods give same results", {
  tis2 <- suppressWarnings(tis(-LLmat, r_eff = r_eff_arr))
  expect_identical(tis1, tis2)

  tisvec <- suppressWarnings(tis(-LLvec, r_eff = r_eff_vec))
  tismat <- suppressWarnings(tis(-LLmat[, 1], r_eff = r_eff_vec))
  expect_identical(tisvec, tismat)

  is2 <- suppressWarnings(sis(-LLmat, r_eff = r_eff_arr))
  expect_identical(is1, is2)

  isvec <- suppressWarnings(sis(-LLvec, r_eff = r_eff_vec))
  ismat <- suppressWarnings(sis(-LLmat[, 1], r_eff = r_eff_vec))
  expect_identical(isvec, ismat)
})



test_that("tis throws correct errors and warnings", {
  # r_eff default no warnings
  expect_silent(tis(-LLarr))
  expect_silent(tis(-LLmat))
  expect_silent(tis(-LLmat[, 1]))

  # r_eff=NULL no warnings
  expect_silent(tis(-LLarr, r_eff = NULL))
  expect_silent(tis(-LLmat, r_eff = NULL))
  expect_silent(tis(-LLmat[,1], r_eff = NULL))

  # r_eff=NA no warnings
  expect_silent(tis(-LLarr, r_eff = NA))
  expect_silent(tis(-LLmat, r_eff = NA))
  expect_silent(tis(-LLmat[,1], r_eff = NA))

  # r_eff default and r_eff=NA give same answer
  expect_equal(
    suppressWarnings(tis(-LLarr)),
                     tis(-LLarr, r_eff = NA)
  )

  # r_eff=NULL and r_eff=NA give same answer
  expect_equal(
    suppressWarnings(tis(-LLarr, r_eff = NULL)),
                     tis(-LLarr, r_eff = NA)
  )

  # r_eff scalar is fine
  expect_silent(tis(-LLarr, r_eff = r_eff_arr[1]))

  # r_eff wrong length is error
  expect_error(tis(-LLarr, r_eff = r_eff_arr[-1]), "one value per observation")

  # r_eff has some NA values causes error
  r_eff_arr[2] <- NA
  expect_error(tis(-LLarr, r_eff = r_eff_arr), "mix NA and not NA values")

  # no NAs or non-finite values allowed
  LLmat[1,1] <- NA
  expect_error(tis(-LLmat), "NAs not allowed in input")

  LLmat[1,1] <- 1
  LLmat[10, 2] <- Inf
  expect_error(tis(-LLmat), "All input values must be finite")

  # no lists allowed
  expect_error(expect_warning(tis(as.list(-LLvec)), "List not allowed as input"))

  # if array, must be 3-D array
  dim(LLarr) <- c(2, 250, 2, 32)
  expect_error(
    tis(-LLarr),
    "length(dim(log_ratios)) == 3 is not TRUE",
    fixed = TRUE
  )
})


test_that("explict test of values for 'sis' and 'tis'", {
  lw <- 1:16
  expect_silent(tis_true <- tis(log_ratios = lw, r_eff = NA))
  expect_equal(as.vector(weights(tis_true, log = TRUE, normalize = FALSE)),
               c(-14.0723, -13.0723, -12.0723, -11.0723, -10.0723, -9.0723, -8.0723, -7.0723, -6.0723, -5.0723, -4.0723, -3.0723, -2.0723, -1.0723, -0.0723, 0.) + 15.07238, tol = 0.001)
  expect_silent(is_true <- sis(log_ratios = lw, r_eff = NA))
  expect_equal(as.vector(weights(is_true, log = TRUE, normalize = FALSE)),
               lw, tol = 0.00001)

  lw <- c(0.7609420, 1.3894140, 0.4158346, 2.5307927, 4.3379119, 2.4159240, 2.2462172, 0.8057697, 0.9333107, 1.5599302)

  expect_silent(tis_true <- tis(log_ratios = lw, r_eff = NA))
  expect_equal(as.vector(weights(tis_true, log = TRUE, normalize = FALSE)),
               c(-2.931, -2.303, -3.276, -1.161,  0, -1.276, -1.446, -2.886, -2.759, -2.132) + 3.692668,
               tol = 0.001)
  expect_silent(is_true <- sis(log_ratios = lw, r_eff = NA))
  expect_equal(as.vector(weights(is_true, log = TRUE, normalize = FALSE)),
               lw, tol = 0.00001)
})


test_that("tis_loo and sis_loo are returned", {
  LLmat <- example_loglik_matrix()
  loo_psis <- suppressWarnings(loo(LLmat, r_eff = NA, is_method = "psis"))
  loo_tis <- suppressWarnings(loo(LLmat, r_eff = NA, is_method = "tis"))
  loo_sis <- suppressWarnings(loo(LLmat, r_eff = NA, is_method = "sis"))

  expect_s3_class(loo_tis, "tis_loo")
  expect_s3_class(loo_sis, "sis_loo")
  expect_s3_class(loo_tis, "importance_sampling_loo")
  expect_s3_class(loo_sis, "importance_sampling_loo")

  expect_output(print(loo_tis), regexp = "tis_loo")
  expect_output(print(loo_sis), regexp = "sis_loo")
})
