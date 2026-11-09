LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()

test_that("example_loglik_array and example_loglik_matrix dimensions ok", {
  dim_arr <- dim(LLarr)
  dim_mat <- dim(LLmat)
  expect_equal(dim_mat[1], dim_arr[1] * dim_arr[2])
  expect_equal(dim_mat[2], dim_arr[3])
})

test_that("example_loglik_array and example_loglik_matrix contain same values", {
  expect_equal(LLmat[1:500, ], LLarr[, 1, ])
  expect_equal(LLmat[501:1000, ], LLarr[, 2, ])
})

test_that("reshaping functions result in correct dimensions", {
  LLmat2 <- llarray_to_matrix(LLarr)
  expect_identical(LLmat2, LLmat)

  LLarr2 <- llmatrix_to_array(LLmat2, chain_id = rep(1:2, each = 500))
  expect_identical(LLarr2, LLarr)
})

test_that("reshaping functions throw correct errors", {
  expect_error(
    llmatrix_to_array(LLmat, chain_id = rep(1:2, times = c(400, 600))),
    regexp = "Not all chains have same number of iterations",
    fixed = TRUE
  )
  expect_error(
    llmatrix_to_array(LLmat, chain_id = rep(1:2, each = 400)),
    regexp = "Number of rows in matrix not equal to length(chain_id)",
    fixed = TRUE
  )
  expect_error(
    llmatrix_to_array(LLmat, chain_id = rep(2:3, each = 500)),
    regexp = "max(chain_id) not equal to the number of chains",
    fixed = TRUE
  )
  expect_error(
    llmatrix_to_array(LLmat, chain_id = rnorm(1000)),
    regexp = "all(chain_id == as.integer(chain_id)) is not TRUE",
    fixed = TRUE
  )
})

test_that("colLogMeanExps(x) = log(colMeans(exp(x))) ", {
  expect_equal(colLogMeanExps(LLmat), log(colMeans(exp(LLmat))))
})

test_that("log-mean-exp helpers handle negative infinity", {
  x <- c(-Inf, log(2), log(3))
  expect_equal(logMeanExp(x), log(5 / 3))
  expect_equal(logMeanExp(rep(-Inf, 3)), -Inf)

  x <- cbind(x, rep(-Inf, 3))
  expect_equal(unname(colLogMeanExps(x)), c(log(5 / 3), -Inf))
})

test_that("validating log-lik objects and functions works", {
  f_ok <- function(data_i, draws) return(NULL)
  f_bad1 <- function(data_i) return(NULL)
  f_bad2 <- function(data, draws) return(NULL)
  expect_equal(validate_llfun(f_ok), f_ok)

  bad_msg <- "Log-likelihood function must have at least the arguments 'data_i' and 'draws'"
  expect_error(validate_llfun(f_bad1), bad_msg)
  expect_error(validate_llfun(f_bad2), bad_msg)
})

test_that("nlist works", {
  a <- 1
  b <- 2
  c <- 3
  nlist_val <- list(nlist(a, b, c), nlist(a, b, c = "tornado"))
  nlist_ans <- list(
    list(a = 1, b = 2, c = 3),
    list(a = 1, b = 2, c = "tornado")
  )
  expect_equal(nlist_val, nlist_ans)
  expect_equal(nlist(a = 1, b = 2, c = 3), list(a = 1, b = 2, c = 3))
})

test_that("loo_cores works", {
  expect_equal(loo_cores(10), 10)
  options(mc.cores = 2)
  expect_equal(loo_cores(getOption("mc.cores", 1)), 2)
  options(mc.cores = 1)

  options(loo.cores = 2)
  expect_warning(expect_equal(loo_cores(10), 2), "deprecated")
  options(loo.cores = NULL)
})

test_that(".collect_warnings() signals each posterior class once", {
  warn_classed <- function(class) {
    warning(structure(
      class = c(class, "posterior_warning", "warning", "condition"),
      list(message = class, call = NULL)
    ))
  }
  caught <- character()
  out <- withCallingHandlers(
    .collect_warnings({
      warn_classed("a")
      warn_classed("a")
      .collect_warnings(warn_classed("b"))
      warn_classed("b")
      warning("other")
      1
    }),
    warning = function(cnd) {
      caught <<- c(caught, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(out, 1)
  expect_equal(caught, c("other", "a (2 times)", "b (2 times)"))
})

test_that(".collect_warnings() keeps the first warning and signals it after an error", {
  warn_a <- function(msg) {
    warning(structure(
      class = c("a", "posterior_warning", "warning", "condition"),
      list(message = msg, call = NULL)
    ))
  }
  w <- capture_warnings(try(.collect_warnings({
    warn_a("first")
    warn_a("second")
    stop("boom")
  }), silent = TRUE))
  expect_equal(w, "first (2 times)")
})

test_that("relative_eff() signals each posterior warning class once", {
  local_mocked_bindings(
    ess_mean = function(x, ...) {
      warning(structure(
        class = c("posterior_warning_x", "posterior_warning", "warning", "condition"),
        list(message = "mocked", call = NULL)
      ))
      1
    },
    .package = "posterior"
  )
  w <- capture_warnings(relative_eff(exp(example_loglik_array()), cores = 1))
  expect_equal(w, "mocked (32 times)")
})

test_that("relative_eff() signals the posterior ESS warning once", {
  skip_if_not(packageVersion("posterior") >= "1.7.1")
  set.seed(1)
  x <- array(rep(c(1, 2), 300) + rnorm(600, 0, 0.1), dim = c(100, 2, 3))
  w <- capture_warnings(relative_eff(x, cores = 1))
  expect_equal(w, "The ESS has been capped to avoid unstable estimates. (3 times)")
})
