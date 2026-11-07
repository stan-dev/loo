# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

# make log lik available to tests
touchstone::pin_assets("touchstone/wine.rds")

# The base branch takes `cores`; the branches with `parallel` run on mirai
# daemons set up beforehand. `sequential` and `two_workers` are the arguments
# that give the same amount of parallelism on either.
touchstone::benchmark_run(
  expr_before_benchmark = {
    suppressPackageStartupMessages(library(loo))
    wine_log_lik_matrix <- readRDS(touchstone::path_pinned_asset(
      "touchstone/wine.rds"
    ))
    matrix_r_eff <- rep(1, ncol(wine_log_lik_matrix))
    uses_daemons <- "parallel" %in% names(formals(loo:::loo.matrix))
    sequential <- if (uses_daemons) list(parallel = FALSE) else list(cores = 1)
    two_workers <- if (uses_daemons) list(parallel = TRUE) else list(cores = 2)
    if (uses_daemons) {
      mirai::daemons(2)
      invisible(mirai::everywhere(library(loo)))
    }
  },
  loo_matrix = {
    suppressWarnings(
      do.call(
        loo,
        c(list(wine_log_lik_matrix, r_eff = matrix_r_eff), sequential)
      )
    )
  },
  loo_matrix_two_workers = {
    suppressWarnings(
      do.call(
        loo,
        c(list(wine_log_lik_matrix, r_eff = matrix_r_eff), two_workers)
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
    uses_daemons <- "parallel" %in% names(formals(loo:::loo.matrix))
    sequential <- if (uses_daemons) list(parallel = FALSE) else list(cores = 1)
    two_workers <- if (uses_daemons) list(parallel = TRUE) else list(cores = 2)
    if (uses_daemons) {
      mirai::daemons(2)
      invisible(mirai::everywhere(library(loo)))
    }
  },
  loo_function = {
    suppressWarnings(
      do.call(
        loo,
        c(
          list(
            wine_llfun,
            data = wine_data,
            draws = wine_log_lik_matrix,
            r_eff = function_r_eff
          ),
          sequential
        )
      )
    )
  },
  loo_function_two_workers = {
    suppressWarnings(
      do.call(
        loo,
        c(
          list(
            wine_llfun,
            data = wine_data,
            draws = wine_log_lik_matrix,
            r_eff = function_r_eff
          ),
          two_workers
        )
      )
    )
  },
  n = 60
)

# create artifacts used downstream in the GitHub Action
touchstone::benchmark_analyze()
