library(loo)
options(mc.cores=1)
options(loo.cores=NULL)
set.seed(123)

context("psis")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = dim(LLarr)[1])
r_eff_arr <- relative_eff(exp(LLarr))
r_eff_vec <- relative_eff(exp(LLvec), chain_id = chain_id)
psis1 <- psis(log_ratios = -LLarr, r_eff = r_eff_arr)

test_that("psis results haven't changed", {
  expect_equal_to_reference(psis1, "reference-results/psis.rds")
})

test_that("psis returns object with correct structure", {
  expect_true(is.psis(psis1))
  expect_false(is.loo(psis1))
  expect_false(is.psis_loo(psis1))

  expect_named(psis1, c("log_weights", "diagnostics"))
  expect_named(psis1$diagnostics, c("pareto_k", "n_eff", "r_eff"))
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
  # r_eff default no warnings
  expect_no_warning(psis(-LLarr))
  expect_no_warning(psis(-LLmat))
  expect_no_warning(psis(-LLmat[, 1]))

  # r_eff=NULL no warnings
  expect_silent(psis(-LLarr, r_eff = NULL))
  expect_silent(psis(-LLmat, r_eff = NULL))
  expect_silent(psis(-LLmat[,1], r_eff = NULL))

  # r_eff=NA disables warnings
  expect_silent(psis(-LLarr, r_eff = NA))
  expect_silent(psis(-LLmat, r_eff = NA))
  expect_silent(psis(-LLmat[,1], r_eff = NA))

  # r_eff default and r_eff=NA give same answer
  expect_equal(
    suppressWarnings(psis(-LLarr)),
    psis(-LLarr, r_eff = NA)
  )

  # r_eff=NULL and r_eff=NA give same answer
  expect_equal(
    suppressWarnings(psis(-LLarr, r_eff=NULL)),
    psis(-LLarr, r_eff = NA)
  )

  # r_eff scalar is fine
  expect_silent(psis(-LLarr, r_eff = r_eff_arr[1]))

  # r_eff non-scalar wrong length is error
  expect_error(psis(-LLarr, r_eff = r_eff_arr[-1]), "one value per observation")

  # r_eff has some NA values causes error
  r_eff_arr[2] <- NA
  expect_error(psis(-LLarr, r_eff = r_eff_arr), "mix NA and not NA values")

  # tail length warnings
  expect_warning(
    psis(-LLarr[1:5,, ]),
    "Not enough tail samples to fit the generalized Pareto distribution"
  )

  # no NAs or non-finite values allowed
  LLmat[1,1] <- NA
  expect_error(psis(-LLmat), "NAs not allowed in input")

  LLmat[1,1] <- 1
  LLmat[10, 2] <- -Inf
  expect_error(psis(-LLmat), "All input values must be finite or -Inf")
  # log ratio of -Inf is allowed
  LLmat[10, 2] <- Inf
  expect_no_error(psis(-LLmat))

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
  expect_no_warning(psis_n_eff.default(w[, 1]))
  expect_no_warning(psis_n_eff.matrix(w))
})


test_that("do_psis_i throws warning if all tail values the same", {
  xx <- c(1,2,3,4,4,4,4,4,4,4,4)
  val <- expect_warning(do_psis_i(xx, tail_len_i = 6), "all tail values are the same")
  expect_equal(val$pareto_k, Inf)
})

test_that("psis_smooth_tail returns original tail values if k is infinite", {
  # skip on M1 Mac until we figure out why this test fails only on M1 Mac
  skip_if(Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64")

  xx <- c(1,2,3,4,4,4,4,4,4,4,4)
  val <- suppressWarnings(psis_smooth_tail(xx, 3))
  expect_equal(val$tail, xx)
  expect_equal(val$k, Inf)
})

