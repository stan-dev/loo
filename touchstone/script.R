# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

# make log lik available to tests
touchstone::pin_assets("touchstone/wine.rds")

# The base branch takes `cores`; branches with `parallel` take `parallel` and
# run on mirai daemons started beforehand. `sequential` and `two_workers` hold
# the arguments that give the same amount of parallelism on either.
setup <- rlang::expr({
  suppressPackageStartupMessages(library(loo))
  wine_log_lik_matrix <- readRDS(touchstone::path_pinned_asset(
    "touchstone/wine.rds"
  ))
  r_eff <- rep(1, ncol(wine_log_lik_matrix))
  wine_data <- data.frame(obs = seq_len(ncol(wine_log_lik_matrix)))
  wine_llfun <- function(data_i, draws) draws[, data_i$obs, drop = FALSE]
  has_parallel <- "parallel" %in% names(formals(loo:::loo.matrix))
  sequential <- if (has_parallel) list(parallel = FALSE) else list(cores = 1)
  two_workers <- if (has_parallel) list(parallel = TRUE) else list(cores = 2)
})

start_daemons <- rlang::expr(
  if (has_parallel) {
    mirai::daemons(2)
    invisible(mirai::everywhere(library(loo)))
  }
)

touchstone::benchmark_run(
  expr_before_benchmark = !!setup,
  loo_matrix = {
    suppressWarnings(
      do.call(loo, c(list(wine_log_lik_matrix, r_eff = r_eff), sequential))
    )
  },
  n = 60
)

touchstone::benchmark_run(
  expr_before_benchmark = {
    !!setup
    !!start_daemons
  },
  loo_matrix_two_workers = {
    suppressWarnings(
      do.call(loo, c(list(wine_log_lik_matrix, r_eff = r_eff), two_workers))
    )
  },
  n = 60
)

touchstone::benchmark_run(
  expr_before_benchmark = !!setup,
  loo_function = {
    suppressWarnings(
      do.call(
        loo,
        c(
          list(
            wine_llfun,
            data = wine_data,
            draws = wine_log_lik_matrix,
            r_eff = r_eff
          ),
          sequential
        )
      )
    )
  },
  n = 60
)

touchstone::benchmark_run(
  expr_before_benchmark = {
    !!setup
    !!start_daemons
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
            r_eff = r_eff
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
