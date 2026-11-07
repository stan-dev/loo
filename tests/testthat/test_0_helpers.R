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


test_that("loo_map re-raises warnings from daemons", {
  expect_warning(
    out <- loo_map(
      1L,
      function(i) {
        warning("from worker")
        i
      },
      parallel = TRUE
    ),
    "from worker"
  )
  expect_equal(out, list(1L))
})

test_that("parallel = TRUE names the missing packages", {
  expect_error(
    check_parallel_packages(c("mirai", "notapackage")),
    "'notapackage' package"
  )
  expect_silent(check_parallel_packages("stats"))
})

test_that("loo_map gives the same result whether arguments are shared or copied", {
  skip_if_not_installed("mori")
  x <- matrix(rnorm(40), 10)
  column_sum <- function(i, x) sum(x[, i])
  copied <- loo_map(1:4, column_sum, x = x, parallel = TRUE, share_bytes = Inf)
  shared <- loo_map(1:4, column_sum, x = x, parallel = TRUE, share_bytes = 0)
  unbatched <- loo_map(1:4, column_sum, x = x, parallel = TRUE, batch = FALSE)
  expect_equal(shared, copied)
  expect_equal(unbatched, copied)
  expect_equal(copied, as.list(colSums(x)), ignore_attr = TRUE)
})

test_that("loo.parallel sets the default of parallel", {
  withr::local_options(loo.parallel = FALSE)
  expect_false(eval(formals(loo.matrix)$parallel))
})
