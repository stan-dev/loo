# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

# These synthetic workloads are large enough to expose real slowdowns in the
# core `loo()` paths, but still short enough to keep PR feedback reasonably fast.
touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))
    matrix_r_eff <- rep(1, ncol(loo:::.wine_log_lik_matrix))
  },
  loo_matrix = {
    suppressWarnings(
      loo(
        loo:::.wine_log_lik_matrix,
        r_eff = matrix_r_eff,
        cores = 1
      )
    )
  },
  n = 10
)

touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))
    function_r_eff <- rep(1, ncol(loo:::.wine_log_lik_matrix))
    wine_data <- data.frame(obs = seq_len(ncol(loo:::.wine_log_lik_matrix)))
    wine_llfun <- function(data_i, draws) draws[, data_i$obs, drop = FALSE]
  },
  loo_function = {
    suppressWarnings(
      loo(
        wine_llfun,
        data = wine_data,
        draws = loo:::.wine_log_lik_matrix,
        r_eff = function_r_eff,
        cores = 1
      )
    )
  },
  n = 10
)

# create artifacts used downstream in the GitHub Action
touchstone::benchmark_analyze()
