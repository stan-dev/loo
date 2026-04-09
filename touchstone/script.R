# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

# These synthetic workloads are large enough to expose real slowdowns in the
# core `loo()` paths, but still short enough to keep PR feedback reasonably fast.
touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))

    matrix_draws <- 2000L
    matrix_obs <- 500L
    n_chains <- 4L
    stopifnot(matrix_draws %% n_chains == 0L)

    set.seed(20260408)
    mu <- stats::rnorm(matrix_draws)
    sigma <- exp(stats::rnorm(matrix_draws, mean = -0.2, sd = 0.15))
    y <- stats::rnorm(matrix_obs, mean = 0.3, sd = 1.2)

    log_lik_matrix <- vapply(
      y,
      FUN = function(y_i) {
        stats::dnorm(y_i, mean = mu, sd = sigma, log = TRUE)
      },
      FUN.VALUE = numeric(matrix_draws)
    )
    matrix_r_eff <- rep(1, matrix_obs)

  },
  loo_matrix = {
    suppressWarnings(loo(log_lik_matrix, r_eff = matrix_r_eff, cores = 1))
  },
  n = 10
)

touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))

    matrix_draws <- 2000L
    matrix_obs <- 500L
    n_chains <- 4L
    stopifnot(matrix_draws %% n_chains == 0L)

    set.seed(20260408)
    mu <- stats::rnorm(matrix_draws)
    sigma <- exp(stats::rnorm(matrix_draws, mean = -0.2, sd = 0.15))
    y <- stats::rnorm(matrix_obs, mean = 0.3, sd = 1.2)

    log_lik_matrix <- vapply(
      y,
      FUN = function(y_i) {
        stats::dnorm(y_i, mean = mu, sd = sigma, log = TRUE)
      },
      FUN.VALUE = numeric(matrix_draws)
    )
    matrix_r_eff <- rep(1, matrix_obs)

    n_iter <- matrix_draws / n_chains
    log_lik_array <- array(NA_real_, dim = c(n_iter, n_chains, matrix_obs))
    for (chain in seq_len(n_chains)) {
      rows <- ((chain - 1L) * n_iter + 1L):(chain * n_iter)
      log_lik_array[, chain, ] <- log_lik_matrix[rows, ]
    }
  },
  loo_array = {
    suppressWarnings(loo(log_lik_array, r_eff = matrix_r_eff, cores = 1))
  },
  n = 10
)

touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))

    matrix_draws <- 2000L
    function_obs <- 250L

    set.seed(20260408)
    mu <- stats::rnorm(matrix_draws)
    sigma <- exp(stats::rnorm(matrix_draws, mean = -0.2, sd = 0.15))
    y <- stats::rnorm(function_obs, mean = 0.3, sd = 1.2)

    function_data <- data.frame(y = y)
    function_draws <- cbind(mu = mu, sigma = sigma)
    function_r_eff <- rep(1, function_obs)
    llfun <- function(data_i, draws) {
      stats::dnorm(
        data_i$y,
        mean = draws[, "mu"],
        sd = draws[, "sigma"],
        log = TRUE
      )
    }
  },
  loo_function = {
    suppressWarnings(
      loo(
        llfun,
        data = function_data,
        draws = function_draws,
        r_eff = function_r_eff,
        cores = 1
      )
    )
  },
  n = 10
)

# create artifacts used downstream in the GitHub Action
touchstone::benchmark_analyze()
