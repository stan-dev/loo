# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

# make log lik available to tests
touchstone::pin_assets("touchstone/wine.rds")

touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))
    wine_log_lik_matrix <- readRDS(touchstone::path_pinned_asset(
      "touchstone/wine.rds"
    ))
    matrix_r_eff <- rep(1, ncol(wine_log_lik_matrix))
  },
  loo_matrix = {
    suppressWarnings(
      loo(
        wine_log_lik_matrix,
        r_eff = matrix_r_eff,
        cores = 1
      )
    )
  },
  n = 60
)

touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))
    wine_log_lik_matrix <- readRDS(touchstone::path_pinned_asset(
      "touchstone/wine.rds"
    ))
    function_r_eff <- rep(1, ncol(wine_log_lik_matrix))
    wine_data <- data.frame(obs = seq_len(ncol(wine_log_lik_matrix)))
    wine_llfun <- function(data_i, draws) draws[, data_i$obs, drop = FALSE]
  },
  loo_function = {
    suppressWarnings(
      loo(
        wine_llfun,
        data = wine_data,
        draws = wine_log_lik_matrix,
        r_eff = function_r_eff,
        cores = 1
      )
    )
  },
  n = 60
)

# create artifacts used downstream in the GitHub Action
touchstone::benchmark_analyze()
