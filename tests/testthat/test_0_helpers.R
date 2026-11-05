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

local_mizu_pool <- function(n = 2L, env = parent.frame()) {
  pool <- mizu::mizu_pool(n)
  previous <- mizu::mizu_set_default_pool(pool)
  withr::defer(
    {
      mizu::mizu_set_default_pool(previous)
      mizu::mizu_pool_stop(pool)
    },
    envir = env
  )
  pool
}

test_that("loo_map on mizu matches lapply for every batch and size", {
  skip_if_not_installed("mizu")
  local_mizu_pool()
  scale_by <- c(10, 20)
  f <- function(i, scale_by) i * scale_by[2]
  expect_true(mizu_pool_is_set())
  for (n in c(1L, 2L, 3L, 7L, 50L, 1001L)) {
    for (batch in c(TRUE, FALSE)) {
      expect_identical(
        loo_map(seq_len(n), f, scale_by = scale_by, parallel = TRUE, batch = batch),
        lapply(seq_len(n), f, scale_by = scale_by)
      )
    }
  }
  expect_identical(
    loo_map(integer(0), f, scale_by = scale_by, parallel = TRUE),
    list()
  )
})

test_that("loo_map on mizu preserves order of unevenly slow elements", {
  skip_if_not_installed("mizu")
  local_mizu_pool()
  slow_early <- function(i) {
    if (i <= 2L) Sys.sleep(0.2)
    i
  }
  expect_identical(
    loo_map(1:20, slow_early, parallel = TRUE, batch = FALSE),
    as.list(1:20)
  )
})

test_that("loo_map on mizu passes mori-shared and plain arguments", {
  skip_if_not_installed("mizu")
  local_mizu_pool()
  big <- matrix(rnorm(2e4), 100)
  lst <- list(a = 1:3, b = "x")
  f <- function(i, big, lst, k) c(sum(big[, i]), length(lst$a), k)
  expect_identical(
    loo_map(1:20, f, big = big, lst = lst, k = 5L, parallel = TRUE),
    lapply(1:20, f, big = big, lst = lst, k = 5L)
  )
})

test_that("loo_map on mizu re-raises every worker warning exactly once", {
  skip_if_not_installed("mizu")
  local_mizu_pool()
  seen <- character()
  out <- withCallingHandlers(
    loo_map(1:6, function(i) {
      warning("w", i)
      i
    }, parallel = TRUE),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(out, as.list(1:6))
  expect_setequal(seen, paste0("w", 1:6))
  expect_length(seen, 6L)
})

test_that("loo_map on mizu stops on the first worker error", {
  skip_if_not_installed("mizu")
  local_mizu_pool()
  expect_error(
    loo_map(1:10, function(i) if (i == 4L) stop("mizu boom") else i, parallel = TRUE),
    "mizu boom"
  )
  expect_error(
    loo_map(1:10, function(i) stop("mizu boom"), parallel = TRUE, batch = FALSE),
    "mizu boom"
  )
  expect_identical(loo_map(1:3, identity, parallel = TRUE), as.list(1:3))
})

test_that("loo_map prefers a mizu pool over mirai daemons and falls back to mirai", {
  skip_if(identical(Sys.getenv("LOO_TEST_BACKEND"), "mizu"))
  skip_if_not_installed("mizu")
  mirai_pids <- unlist(loo_map(1:40, function(i) Sys.getpid(), parallel = TRUE, batch = FALSE))
  expect_false(mizu_pool_is_set())

  pool <- local_mizu_pool()
  mizu_pids <- unlist(loo_map(1:40, function(i) Sys.getpid(), parallel = TRUE, batch = FALSE))
  expect_true(mizu_pool_is_set())
  expect_length(intersect(mirai_pids, mizu_pids), 0L)
  expect_false(Sys.getpid() %in% mizu_pids)

  mizu::mizu_set_default_pool(NULL)
  expect_false(mizu_pool_is_set())
  again <- unlist(loo_map(1:40, function(i) Sys.getpid(), parallel = TRUE, batch = FALSE))
  expect_true(all(again %in% mirai_pids))
  mizu::mizu_set_default_pool(pool)
})

test_that("parallel = FALSE ignores a mizu pool", {
  skip_if_not_installed("mizu")
  local_mizu_pool()
  pids <- unlist(loo_map(1:3, function(i) Sys.getpid(), parallel = FALSE))
  expect_true(all(pids == Sys.getpid()))
})

test_that("psis, relative_eff and loo give identical results on mizu, mirai and sequentially", {
  skip_if_not_installed("mizu")
  set.seed(1)
  S <- 400
  N <- 40
  ll <- matrix(rnorm(S * N, -1, 0.5), S, N)
  ll[, 1:3] <- ll[, 1:3] * c(1, 4)
  arr <- array(ll, c(100, 4, N))
  llfun <- function(data_i, draws) dnorm(data_i$y, draws[, 1], 1, log = TRUE)
  dat <- data.frame(y = rnorm(N))
  dr <- cbind(rnorm(S))

  compute <- function(parallel) {
    suppressWarnings(list(
      r_eff = relative_eff(exp(arr), parallel = parallel),
      psis = psis(-ll, r_eff = rep(1, N), parallel = parallel),
      tis = tis(-ll, r_eff = rep(1, N), parallel = parallel),
      loo_array = loo(arr, r_eff = rep(1, N), parallel = parallel),
      loo_fun = loo(llfun, data = dat, draws = dr, r_eff = rep(1, N), parallel = parallel),
      r_eff_fun = relative_eff(
        function(data_i, draws, ...) exp(llfun(data_i, draws)),
        chain_id = rep(1:2, each = S / 2), data = dat, draws = dr, parallel = parallel
      )
    ))
  }
  sequential <- compute(FALSE)
  on_mirai <- compute(TRUE)
  local_mizu_pool()
  on_mizu <- compute(TRUE)
  expect_identical(on_mirai, sequential)
  expect_identical(on_mizu, sequential)
})
