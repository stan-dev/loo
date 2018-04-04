library(loo)
suppressPackageStartupMessages(library(rstanarm))
options(mc.cores=1)
options(loo.cores=NULL)
set.seed(123)

context("psis")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = dim(LLarr)[1])
r_eff_arr <- relative_eff(exp(-LLarr))
r_eff_vec <- relative_eff(exp(-LLvec), chain_id = chain_id)
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

test_that("psis throws correct errors and warnings", {
  # r_eff warnings
  expect_warning(psis(-LLarr), "Relative effective sample sizes")
  expect_warning(psis(-LLmat), "Relative effective sample sizes")
  expect_warning(psis(-LLmat[, 1]), "Relative effective sample sizes")

  # tail length warnings
  expect_warning(
    psis(-LLarr[1:5,, ]),
    "Not enough tail samples to fit the generalized Pareto distribution"
  )

  # no NAs or non-finite values allowed
  LLmat[1,1] <- NA
  expect_error(psis(-LLmat), "NAs not allowed in input")

  LLmat[1,1] <- 1
  LLmat[10, 2] <- Inf
  expect_error(psis(-LLmat), "All input values must be finite")

  # no lists allowed
  expect_error(expect_warning(psis(as.list(-LLvec))), "List not allowed as input")

  # if array, must be 3-D array
  dim(LLarr) <- c(2, 250, 2, 32)
  expect_error(
    psis(-LLarr),
    "length(dim(log_ratios)) == 3 is not TRUE",
    fixed = TRUE
  )
})

test_that("throw_tail_length_warnings gives correct output", {
  expect_silent(throw_tail_length_warnings(10))
  expect_equal(throw_tail_length_warnings(10), 10)
  expect_warning(throw_tail_length_warnings(1), "Not enough tail samples")
  expect_warning(throw_tail_length_warnings(c(1, 10, 2)),
                 "Skipping the following columns: 1, 3")
  expect_warning(throw_tail_length_warnings(rep(1, 21)),
                 "11 more not printed")
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

test_that("relative_eff methods works properly", {
  expect_equal(relative_eff.default(exp(LLmat[, 1]), chain_id),
               mcmc_n_eff(exp(LLarr[, , 1])) / 1000)
  expect_equal(relative_eff.matrix(exp(LLmat), chain_id),
               apply(exp(LLarr), 3, mcmc_n_eff) / 1000)
  expect_equal(relative_eff.array(exp(LLarr)),
               apply(exp(LLarr), 3, mcmc_n_eff) / 1000)
  expect_equal(relative_eff.array(exp(LLarr)),
               apply(exp(LLarr), 3, mcmc_n_eff) / 1000)

  expect_equal(relative_eff(exp(LLarr)), relative_eff(exp(LLarr), cores = 2))
})

test_that("relative_eff function method works properly", {
  dat <- data.frame(y = mtcars$mpg, x = mtcars$wt)
  utils::capture.output(
    fit <- suppressWarnings(
      rstanarm::stan_glm(y ~ x, data = dat, chains = 2, iter = 250, seed = 123)
    )
  )
  draws <- as.matrix(fit)
  chain_id <- rep(1:2, each = 125)
  likfun <- function(data_i, draws) {
    dnorm(data_i$y, draws[, 1] + draws[, 2] * data_i$x, draws[, 3])
  }

  r_eff_fn <- relative_eff.function(likfun, chain_id = chain_id,
                                    data = dat, draws = draws)
  r_eff_mat <- relative_eff.matrix(exp(log_lik(fit)), chain_id)
  expect_equal(r_eff_fn, r_eff_mat)

  r_eff_fn_cores <- relative_eff.function(likfun, chain_id = chain_id,
                                          data = dat, draws = draws, cores = 2)
  expect_equal(r_eff_fn, r_eff_fn_cores)
})
