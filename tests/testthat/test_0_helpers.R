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

test_that("loo_map runs on a default mizu pool and matches lapply", {
  skip_if_not_installed("mizu")
  pool <- mizu::mizu_pool(2L)
  withr::defer(mizu::mizu_pool_stop(pool))
  previous <- mizu::mizu_set_default_pool(pool)
  withr::defer(mizu::mizu_set_default_pool(previous))

  scale_by <- c(10, 20)
  f <- function(i, scale_by) i * scale_by[2]
  expect_true(mizu_pool_is_set())
  expect_equal(
    loo_map(1:50, f, scale_by = scale_by, parallel = TRUE),
    lapply(1:50, f, scale_by = scale_by)
  )
  expect_equal(
    loo_map(1:7, f, scale_by = scale_by, parallel = TRUE, batch = FALSE),
    lapply(1:7, f, scale_by = scale_by)
  )
  expect_equal(loo_map(integer(0), f, scale_by = scale_by, parallel = TRUE), list())
})

test_that("loo_map re-raises warnings and errors from mizu workers", {
  skip_if_not_installed("mizu")
  pool <- mizu::mizu_pool(2L)
  withr::defer(mizu::mizu_pool_stop(pool))
  previous <- mizu::mizu_set_default_pool(pool)
  withr::defer(mizu::mizu_set_default_pool(previous))

  expect_warning(
    out <- loo_map(1L, function(i) {
      warning("from mizu worker")
      i
    }, parallel = TRUE),
    "from mizu worker"
  )
  expect_equal(out, list(1L))
  expect_error(
    loo_map(1:3, function(i) stop("mizu boom"), parallel = TRUE),
    "mizu boom"
  )
})

test_that("loo gives the same result on mizu, mirai and sequentially", {
  skip_if_not_installed("mizu")
  set.seed(1)
  ll <- matrix(rnorm(400 * 30, -1, 0.5), 400, 30)
  sequential <- suppressWarnings(loo(ll, r_eff = rep(1, 30), parallel = FALSE))
  on_mirai <- suppressWarnings(loo(ll, r_eff = rep(1, 30), parallel = TRUE))

  pool <- mizu::mizu_pool(2L)
  withr::defer(mizu::mizu_pool_stop(pool))
  previous <- mizu::mizu_set_default_pool(pool)
  withr::defer(mizu::mizu_set_default_pool(previous))
  on_mizu <- suppressWarnings(loo(ll, r_eff = rep(1, 30), parallel = TRUE))

  expect_equal(on_mirai$estimates, sequential$estimates)
  expect_equal(on_mizu$estimates, sequential$estimates)
  expect_equal(on_mizu$diagnostics, sequential$diagnostics)
})
