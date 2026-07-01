options(mc.cores = 1)
set.seed(123)

# Make sure no daemon pool leaks in from another test file.
mirai::daemons(0)

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
chain_id <- rep(1:2, each = dim(LLarr)[1])
r_eff <- relative_eff(exp(LLarr))

# Shared data for the function-method end-to-end checks.
set.seed(1)
S_fn <- 200
N_fn <- 30
draws_fn <- cbind(mu = rnorm(S_fn), sigma = abs(rnorm(S_fn)) + 0.5)
data_fn <- data.frame(y = rnorm(N_fn))
llfun_test <- function(data_i, draws, ...) {
  dnorm(data_i$y, mean = draws[, "mu"], sd = draws[, "sigma"], log = TRUE)
}


# Pool-introspection helpers -------------------------------------------------

test_that("loo_has_pool() and loo_pool_is_local() detect a local pool", {
  mirai::daemons(0)
  expect_false(loo:::loo_has_pool())
  expect_false(loo:::loo_pool_is_local())

  skip_on_cran()
  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  expect_true(loo:::loo_has_pool())
  expect_true(loo:::loo_pool_is_local())
  # Chunking uses the connected daemon count, not the requested cores.
  expect_equal(loo:::loo_n_workers(1), 2L)
})

test_that("loo_pool_is_local() is FALSE for a tcp pool (remote-safety gate)", {
  skip_on_cran()
  mirai::daemons(0)
  mirai::daemons(n = 2, url = mirai::local_url(tcp = TRUE))
  on.exit(mirai::daemons(0), add = TRUE)
  # The locality gate reads the configured transport URL (available
  # immediately, regardless of connection timing). tcp:// may be a remote/SSH
  # pool, so shared memory must not be assumed.
  expect_false(loo:::loo_pool_is_local())
})


# loo_map() ------------------------------------------------------------------

test_that("loo_map() runs serially when no pool is available", {
  mirai::daemons(0)
  res <- loo:::loo_map(1:5, function(x, m) x * m, m = 2, cores = 4)
  expect_identical(res, as.list((1:5) * 2))
})

test_that("loo_map() uses a connected pool regardless of cores", {
  skip_on_cran()
  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  # cores = 1 no longer forces serial: a connected pool always wins.
  res <- loo:::loo_map(1:5, function(x, m) x * m, m = 3, cores = 1)
  expect_identical(res, as.list((1:5) * 3))
})

test_that("loo_map() parallel matches serial and preserves order", {
  skip_on_cran()
  worker <- function(i, mat, add) sum(mat[, i]) + add
  mat <- matrix(as.numeric(1:60), nrow = 6) # 6 x 10
  N <- ncol(mat)
  expected <- lapply(seq_len(N), worker, mat = mat, add = 100)

  mirai::daemons(3)
  on.exit(mirai::daemons(0), add = TRUE)

  # broadcast object shared via mori on a local pool; both chunk strategies
  res_auto <- loo:::loo_map(
    seq_len(N), worker, add = 100, cores = 3,
    broadcast = list(mat = mat), chunk = "auto"
  )
  res_never <- loo:::loo_map(
    seq_len(N), worker, add = 100, cores = 3,
    broadcast = list(mat = mat), chunk = "never"
  )
  expect_identical(res_auto, expected)
  expect_identical(res_never, expected)
})

test_that("loo_map() works when there are more workers than elements", {
  skip_on_cran()
  mirai::daemons(4)
  on.exit(mirai::daemons(0), add = TRUE)
  res <- loo:::loo_map(1:2, function(x) x + 1L, cores = 4)
  expect_identical(res, list(2L, 3L))
})

test_that("loo_map() propagates worker errors", {
  skip_on_cran()
  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  expect_error(
    loo:::loo_map(1:4, function(x) if (x == 3L) stop("boom") else x, cores = 2),
    "boom"
  )
})


# End-to-end: importance sampling --------------------------------------------

test_that("psis() parallel equals serial", {
  skip_on_cran()
  ps_serial <- suppressWarnings(psis(-LLmat, r_eff = r_eff, cores = 1))

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  ps_parallel <- suppressWarnings(psis(-LLmat, r_eff = r_eff, cores = 2))

  expect_equal(ps_serial$log_weights, ps_parallel$log_weights)
  expect_equal(ps_serial$diagnostics, ps_parallel$diagnostics)
})

test_that("tis() and sis() parallel equal serial", {
  skip_on_cran()
  tis_serial <- suppressWarnings(tis(-LLmat, r_eff = r_eff, cores = 1))
  sis_serial <- suppressWarnings(sis(-LLmat, r_eff = r_eff, cores = 1))

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  tis_parallel <- suppressWarnings(tis(-LLmat, r_eff = r_eff, cores = 2))
  sis_parallel <- suppressWarnings(sis(-LLmat, r_eff = r_eff, cores = 2))

  expect_equal(tis_serial$log_weights, tis_parallel$log_weights)
  expect_equal(sis_serial$log_weights, sis_parallel$log_weights)
})


# End-to-end: loo() function method (broadcast draws/data) -------------------

test_that("loo.function parallel equals serial", {
  skip_on_cran()
  loo_serial <- suppressWarnings(
    loo(llfun_test, data = data_fn, draws = draws_fn, cores = 1)
  )

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  loo_parallel <- suppressWarnings(
    loo(llfun_test, data = data_fn, draws = draws_fn, cores = 2)
  )

  expect_equal(loo_serial$pointwise, loo_parallel$pointwise)
  expect_equal(loo_serial$estimates, loo_parallel$estimates)
})

test_that("loo.function reuses an existing (user-configured) pool", {
  skip_on_cran()
  loo_serial <- suppressWarnings(
    loo(llfun_test, data = data_fn, draws = draws_fn, cores = 1)
  )

  # User sets up the pool themselves; loo should reuse it untouched.
  mirai::daemons(3)
  on.exit(mirai::daemons(0), add = TRUE)
  loo_reuse <- suppressWarnings(
    loo(llfun_test, data = data_fn, draws = draws_fn, cores = 2)
  )
  # Pool is still alive after the call (loo did not tear it down).
  expect_true(loo:::loo_has_pool())
  expect_equal(loo_serial$pointwise, loo_reuse$pointwise)
})


# End-to-end: relative_eff ---------------------------------------------------

test_that("relative_eff() array and function methods are parallel-invariant", {
  skip_on_cran()
  re_arr_serial <- relative_eff(exp(LLarr), cores = 1)
  re_fn_serial <- relative_eff(
    llfun_test, chain_id = rep(1, S_fn),
    data = data_fn, draws = draws_fn, cores = 1
  )

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  re_arr_parallel <- relative_eff(exp(LLarr), cores = 2)
  re_fn_parallel <- relative_eff(
    llfun_test, chain_id = rep(1, S_fn),
    data = data_fn, draws = draws_fn, cores = 2
  )

  expect_equal(re_arr_serial, re_arr_parallel)
  expect_equal(re_fn_serial, re_fn_parallel)
})


# End-to-end: loo_subsample --------------------------------------------------

test_that("loo_subsample() parallel equals serial", {
  skip_on_cran()
  # Reset RNG before each call so the same subsample is drawn.
  set.seed(4242)
  ss_serial <- suppressWarnings(loo_subsample(
    llfun_test, data = data_fn, draws = draws_fn,
    observations = 20, loo_approximation = "plpd", cores = 1
  ))

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  set.seed(4242)
  ss_parallel <- suppressWarnings(loo_subsample(
    llfun_test, data = data_fn, draws = draws_fn,
    observations = 20, loo_approximation = "plpd", cores = 2
  ))

  expect_equal(ss_serial$estimates, ss_parallel$estimates)
  expect_equal(ss_serial$pointwise, ss_parallel$pointwise)
})


# End-to-end: loo_model_weights (single pool across K models) ----------------

test_that("loo_model_weights() parallel equals serial", {
  skip_on_cran()
  set.seed(11)
  ll_list <- list(
    matrix(rnorm(200 * 25), nrow = 200),
    matrix(rnorm(200 * 25), nrow = 200),
    matrix(rnorm(200 * 25), nrow = 200)
  )
  wts_serial <- suppressWarnings(
    loo_model_weights(ll_list, method = "stacking", cores = 1)
  )

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  wts_parallel <- suppressWarnings(
    loo_model_weights(ll_list, method = "stacking", cores = 2)
  )

  expect_equal(as.numeric(wts_serial), as.numeric(wts_parallel),
               tolerance = 1e-6)
})

# Persistent session pool (loo.daemons / LOO_DAEMONS) -----------------------

test_that("loo_persist_config() resolves option/env var and rejects bad values", {
  # Start from a known-clean configuration and restore it afterwards.
  old_opt <- options(loo.daemons = NULL)
  on.exit(options(old_opt), add = TRUE)
  old_env <- Sys.getenv("LOO_DAEMONS", unset = NA)
  Sys.unsetenv("LOO_DAEMONS")
  on.exit(
    if (!is.na(old_env)) {
      Sys.setenv(LOO_DAEMONS = old_env)
    } else {
      Sys.unsetenv("LOO_DAEMONS")
    },
    add = TRUE
  )

  # Unset -> feature off.
  expect_identical(loo:::loo_persist_config(), NA_integer_)

  # Environment variable is parsed when the option is unset.
  Sys.setenv(LOO_DAEMONS = "3")
  expect_identical(loo:::loo_persist_config(), 3L)

  # Option takes precedence over the environment variable.
  options(loo.daemons = 4)
  expect_identical(loo:::loo_persist_config(), 4L)

  # 0/1, non-integer, and garbage values all disable the feature (no error).
  options(loo.daemons = 1)
  expect_identical(loo:::loo_persist_config(), NA_integer_)
  options(loo.daemons = 0)
  expect_identical(loo:::loo_persist_config(), NA_integer_)
  options(loo.daemons = 2.5)
  expect_identical(loo:::loo_persist_config(), NA_integer_)
  options(loo.daemons = "garbage")
  expect_identical(suppressWarnings(loo:::loo_persist_config()), NA_integer_)
})

test_that("persistent pool is created lazily and reused across calls", {
  skip_on_cran()
  mirai::daemons(0)
  expect_false(loo:::loo_has_pool())

  old_opt <- options(loo.daemons = 2)
  on.exit(options(old_opt), add = TRUE)
  on.exit(mirai::daemons(0), add = TRUE)

  ps_serial <- suppressWarnings(psis(-LLmat, r_eff = r_eff, cores = 1))
  # cores = 1 work must not spin up the persistent pool.
  expect_false(loo:::loo_has_pool())

  # First parallel call creates the pool and leaves it warm.
  ps1 <- suppressWarnings(psis(-LLmat, r_eff = r_eff, cores = 2))
  expect_true(loo:::loo_has_pool())
  expect_equal(loo:::loo_n_workers(2), 2L)

  # Second call reuses the same pool (still 2 connected daemons).
  ps2 <- suppressWarnings(psis(-LLmat, r_eff = r_eff, cores = 2))
  expect_true(loo:::loo_has_pool())
  expect_equal(loo:::loo_n_workers(2), 2L)

  # Results match the serial computation.
  expect_equal(ps1$log_weights, ps_serial$log_weights)
  expect_equal(ps2$log_weights, ps_serial$log_weights)
})

test_that("a user-configured pool takes precedence over loo.daemons", {
  skip_on_cran()
  mirai::daemons(0)
  old_opt <- options(loo.daemons = 4)
  on.exit(options(old_opt), add = TRUE)

  # User sets up their own pool; loo must reuse it untouched and not replace
  # it with a persistent pool of the configured size.
  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  ps <- suppressWarnings(psis(-LLmat, r_eff = r_eff, cores = 2))
  expect_true(loo:::loo_has_pool())
  expect_equal(loo:::loo_n_workers(2), 2L)
})

test_that("with_loo_daemons() informs (once) that cores is ignored with a pool", {
  skip_on_cran()
  mirai::daemons(0)
  # Reset the once-per-session guard so this test is order-independent.
  loo:::.loo_internal$informed_cores_ignored <- NULL
  on.exit(loo:::.loo_internal$informed_cores_ignored <- NULL, add = TRUE)

  # No pool connected: cores = 1 runs serially without a message.
  expect_silent(loo:::with_loo_daemons(1, 42))

  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # Pool connected and cores = 1: informs once that cores is ignored, and the
  # value still comes back (work runs on the pool).
  expect_message(
    out <- loo:::with_loo_daemons(1, 42),
    "'cores' is ignored"
  )
  expect_identical(out, 42)

  # Message is emitted at most once per session.
  expect_silent(loo:::with_loo_daemons(1, 42))

  # cores > 1 with a pool never triggers the message (parallel was requested).
  loo:::.loo_internal$informed_cores_ignored <- NULL
  expect_silent(loo:::with_loo_daemons(2, 42))
})

test_that("a persistent-pool child process exits cleanly (no orphans)", {
  skip_on_cran()
  skip_on_os("windows")
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if_not(file.exists(rscript))

  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(
    c(
      "library(loo)",
      "LLarr <- example_loglik_array()",
      "r_eff <- relative_eff(exp(LLarr))",
      "invisible(suppressWarnings(loo(LLarr, r_eff = r_eff, cores = 2)))",
      "cat('LOO_CHILD_OK\\n')"
    ),
    script
  )
  out <- suppressWarnings(system2(
    rscript,
    c("--vanilla", shQuote(script)),
    stdout = TRUE,
    stderr = TRUE,
    env = "LOO_DAEMONS=2"
  ))
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(as.integer(status), 0L))
  expect_true(any(grepl("LOO_CHILD_OK", out)))
})

# Final safety net in case any test above exited early with a live pool.
mirai::daemons(0)
