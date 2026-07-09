# Parallelization with mirai in other packages

This document surveys how R packages that depend on **mirai** integrate it into 
their workflows. The goal is to inform design choices for **loo** 
(user-facing API) by comparing established patterns.

## Summary

Inspected packages fall into two broad approaches:

| Approach | User experience | Examples |
|----------|-------------------|----------|
| **Wrapper API** | A `workers` (or similar) argument on a domain function; mirai stays internal | bakerrr, bregr, pliman, GitStats |
| **Direct mirai** | Users call `mirai::daemons()` themselves, then use package functions unchanged | slideimp, worldmet, rush |
| **Abstraction layer** | Another parallel framework (future, crew) uses mirai as a backend | chopin, crew |

Most wrapper-style packages call `mirai::daemons(n)` before parallel work and 
`mirai::daemons(0)` (or `on.exit`) afterward. Several also expose convenience 
helpers so users never touch mirai directly.

For **loo**, two complementary options are worth considering:

1. **Use mirai directly and document workflow** in a vignette

```
mirai::daemons(n)

loo(...)

mirai::daemons(0)
```
2. **Provide a thin convenience wrapper** for users who want parallel 
map-style execution without managing daemons:

```
loo_mirai(fun, args_list, n_daemons = NULL)
```

Third path (currently implemented): 
+ lazy session-scoped pools via `options(loo.daemons)` / `LOO_DAEMONS`, without 
a per-call `workers` argument or a `loo_daemons()` wrapper.


## bakerrr

- **GitHub:** https://github.com/anirbanshaw24/bakerrr
- **Description:** {bakerrr} provides a clean, modern interface for running 
background parallel jobs using S7 classes, mirai daemon(s), and callr process 
management.
- **Pattern:** High-level job API; `n_daemons` controls pool size 
(default: `ceiling(cores/5)`).

### Mirai example usage

- `n_daemons` — number of parallel workers (default: ceiling of cores/5).

#### API

```
# Define your function
compute_sum <- function(x, y) {
  Sys.sleep(1)  # Simulate work
  x + y
}

# Create argument lists for each job
args_list <- list(
  list(x = 1, y = 2),
  list(x = 3, y = 4),
  list(x = 5, y = 6),
  list(x = 7, y = 8)
)

# Create and run bakerrr job
job <- bakerrr::bakerrr(
  fun = compute_sum,
  args_list = args_list,
  n_daemons = 2
) |> 
  bakerrr::run_jobs(wait_for_results = TRUE)
```

#### Backend

```
mirai::daemons(n_daemons)
  results <- purrr::imap(
    jobs,
    purrr::in_parallel(
      \(x, y) {
        tryCatch({
          do.call(x$fun, x$args)
        },
        error = function(e) {
          glue::glue(
            "Error in purrr::in_parallel: {conditionMessage(e)}"
          )
        })
      }
    )
  )
mirai::daemons(0)
```

## bregr

- **GitHub:** https://github.com/WangLabCSU/bregr
- **Description:** The bregr package provides batch regression modeling 
in R, enabling you to run hundreds of models simultaneously with a clean, 
intuitive workflow.
- **Pattern:** `n_workers` on the pipeline function; uses 
`with(daemons(n_workers), mirai_map(...))` for parallel paths and serial `map()`
 when `n_workers == 1`.

### Mirai example usage

- `n_workers` — integer number of workers to launch (default `1L`). When `> 1`, 
modeling runs in parallel in the background.

#### API

```
lung <- survival::lung |>
  dplyr::filter(ph.ecog != 3)
lung$ph.ecog <- factor(lung$ph.ecog)

mds_p <- br_pipeline(
  lung,
  y = c("time", "status"),
  x = colnames(lung)[6:10],
  x2 = c("age", "sex"),
  method = "coxph",
  n_workers = 3
)
#> exponentiate estimates of model(s) constructed from coxph method at default
#> ■■■■■■■                           20% | ETA: 13s
#> 
#> ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s
```

#### Backend

```
if (n_workers > 1) {
    res <- with(
      mirai::daemons(n_workers),
      {
        mp <- mirai::mirai_map(
          ms, runner_,
          .args = list(
            data = data, dots = dots,
            opts = options()
          )
        )
        mp[.progress]
      }
    )
  } else {
    res <- map(
      ms, runner_,
      data = data, dots = dots,
      opts = options()
    )
  }
```

## chopin

- **GitHub:** https://github.com/ropensci/chopin
- **Description:** This package automates parallelization in spatial operations 
with chopin functions as well as sf/terra functions.
- **Pattern:** Offers both **future.mirai** 
(`future::plan(future.mirai::mirai_multisession)`) and a `par_grid_mirai()` 
that expects users to manage `mirai::daemons()` directly. 

### Mirai example usage

#### API

```
# Users always need to register multiple CPU threads (logical cores) for 
# parallelization. 
future::plan(future.mirai::mirai_multisession, workers = 4L)

system.time(
  ncpoints_srtm_mthr <-
    par_grid(
      grids = compregions,
      fun_dist = extract_at,
      x = srtm_ras,
      y = ncpoints,
      id = "pid",
      radius = 1e4L,
      .standalone = FALSE
    )
)

# The same workflow operates on mirai dispatchers.
mirai::daemons(4L)

system.time(
  ncpoints_srtm_mthr <-
    par_grid_mirai(
      grids = compregions,
      fun_dist = extract_at,
      x = srtm_ras,
      y = ncpoints,
      id = "pid",
      radius = 1e4L,
      .standalone = FALSE
    )
)

mirai::daemons(0L)
```

## crew

- **GitHub:** https://github.com/wlandau/crew/
- **Description:** crew is a distributed computing framework with a centralized
 interface and auto-scaling.
- **Pattern:** Controller object (`crew_controller_local`) manages mirai client
 lifecycle.

### Mirai example usage

#### API

```
controller <- crew_controller_local(
  name = "example",
  workers = 2,
  seconds_idle = 10
)

# create the mirai client
controller$start()

# submit a new task
controller$push(name = "get pid", command = ps::ps_pid())

# return completed task
task <- controller$pop()

# clean-up resources
controller$terminate()
```
#### Backend

```
crew_client <- function(...) {
  ...
  #' @description Start listening for workers on the available sockets.
    #' @return `NULL` (invisibly).
    start = function() {
      if (!is.null(.started) && .started) {
        return(invisible())
      }
      existing_url <- mirai::nextget("url", .compute = .profile)
      crew_assert(
        is.null(existing_url),
        message = sprintf(
          paste(
            "{mirai} compute profile '%s' is already active at URL %s.",
            "Each {crew} controller must have a unique compute profile",
            "that does not conflict with",
            "an existing call to mirai::daemons()."
          ),
          .profile,
          existing_url
        )
      )
      mirai::daemons(
        url = .tls$url(host = .host, port = .port),
        dispatcher = TRUE,
        seed = NULL,
        serial = .serialization,
        tls = .tls$client(),
        pass = .tls$password,
        .compute = .profile
      )
      .url <<- mirai::nextget("url", .compute = .profile)
      .relay$set_from(mirai::nextget(x = "cv", .compute = .profile))
      .relay$start()
      .started <<- TRUE
      invisible()
    },
    #' @description Stop the mirai client and disconnect from the
    #'   worker websockets.
    #' @return `NULL` (invisibly).
    terminate = function() {
      mirai::daemons(n = 0L, .compute = .profile)
      .relay$terminate()
      .url <<- NULL
      .started <<- FALSE
      invisible()
    },
    ...
}
```

## GitStats

- **GitHub:** https://github.com/r-world-devs/GitStats/tree/4942e5e569d18ea7064bba0901561e26265c0f8b
- **Description:** Obtain standardized data from multiple 'Git' services, 
including 'GitHub' and 'GitLab'.

### Mirai example usage

#### API

```
git_stats |>
  set_parallel()
  
```
#### Backend

```
set_parallel = function(workers = 10L) {
      if (isFALSE(workers) || identical(workers, 0L) || identical(workers, 0)) {
        status <- mirai::daemons(0)
        private$settings$parallel <- FALSE
        cli::cli_alert_info("Parallel processing disabled.")
      } else {
        if (isTRUE(workers)) {
          workers <- parallel::detectCores(logical = FALSE)
          if (is.na(workers) || workers < 2L) workers <- 2L
        }
        status <- mirai::daemons(workers)
        ns <- asNamespace("GitStats")
        ns_objects <- as.list(ns, all.names = TRUE)
        do.call(mirai::everywhere, c(list(quote({})), ns_objects))
        private$settings$parallel <- TRUE
        cli::cli_alert_success(
          "Parallel processing enabled with {workers} workers."
        )
      }
      return(invisible(status))
    },
  ...
```

## pliman

- **GitHub:** https://github.com/NEPEM-UFSC/pliman
- **Description:** Tools for both single and batch image manipulation and 
analysis and phytopathometry.
- **Pattern:** `parallel` + `workers` arguments on `measure_injury()`; 
`on.exit(daemons(0))` ensures cleanup even on error. Default worker count is 
`trunc(detectCores() * 0.3)`.

### Mirai example usage

#### Backend

```
measure_injury <- function(..., parallel = FALSE, workers = NULL, verbose = TRUE) {
  ...
  init_time <- Sys.time()

  if (parallel) {
    # setup
    # capture start time
    init_time <- Sys.time()

    # decide number of workers
    nworkers <- if (is.null(workers)) trunc(parallel::detectCores() * 0.3) else workers

    # start mirai daemons
    mirai::daemons(nworkers)
    on.exit(mirai::daemons(0), add = TRUE)
    ...
    # run all images in parallel
    raw <- mirai::mirai_map(
      .x = files,
      .f = process_image
    )[.progress]
  }
  ...
}
```

## rush

- **GitHub:** https://github.com/mlr-org/rush
- **Description:** Package designed to solve large-scale problems 
asynchronously across a distributed network.
- **Pattern:** user starts `mirai::daemons()` and `rush$start_workers()` 
separately. 

### Mirai example usage

#### API

```
config = redux::redis_config()
rush = rsh(network_id = "my_network", config = config)

mirai::daemons(n = 2L)

rush$start_workers(
  worker_loop = worker_loop,
  n_workers = 2L)
  
rush$fetch_finished_tasks()
  
rush$stop_workers(type = "kill")
```

## slideimp

- **GitHub:** https://github.com/hhp94/slideimp
- **Description:** KNN and PCA Imputation of Numeric Data.
- **Pattern:** Documentation tells users to call `mirai::daemons(n)` before 
tuning functions and `mirai::daemons(0)` when finished. No `workers` argument 
on the package API — parallelization is entirely user-driven.

### Mirai example usage

#### API

```
mirai::daemons(2) # 2 Cores

# PCA imputation.
pca_params <- data.frame(ncp = c(1, 5))
# For machines with multi-threaded BLAS, turn on `pin_blas = TRUE`
tune_pca <- tune_imp(obj, parameters = pca_params, .f = "pca_imp", n_reps = 10, 
num_na = 50)

mirai::daemons(0) # Close daemons
```

## worldmet

- **GitHub:** https://github.com/openair-project/worldmet
- **Description:** Open source tools to access NOAA meteorological observations.
- **Pattern:** Same as slideimp: set daemons once per session; import functions 
work unchanged when a pool is active.

### Mirai example usage

#### API

```
# set workers - once per session
mirai::daemons(4)

# import lots of data - NB: no change in import_ghcn_hourly()!
big_met <- import_ghcn_hourly(code = "UKI0000EGLL", year = 2010:2025)
```
