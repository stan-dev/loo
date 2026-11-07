args <- commandArgs(TRUE)
impl <- args[1]
book_ll <- args[2]
loo_data <- args[3]
lib <- path.expand(file.path("~/lx", if (impl == "base") "libb" else "libc"))
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = lib)
suppressMessages(library(loo))

read_book <- function(names) lapply(names, function(n) list(ll = readRDS(file.path(book_ll, paste0(n, ".rds")))))
read_loo <- function(file, entries) {
  d <- readRDS(file.path(loo_data, file))
  lapply(entries, function(e) list(ll = d[[e]], chain_id = d$chain_id))
}
milk_kline <- function() {
  d <- readRDS(file.path(loo_data, "real_milk_kline.rds"))
  unname(lapply(d, function(set) lapply(set$ll, function(ll) list(ll = ll, chain_id = set$chain_id))))
}

workloads <- list(
  `bayesian-workflow/roaches` = list(read_book(c("roaches_poisson", "roaches_nb", "roaches_poisson_vi"))),
  `bayesian-workflow/arsenic_loo` = list(read_book(c("arsenic_3", "arsenic_3l"))),
  `loo2-example` = list(read_loo("real_roach.rds", c("ll1", "ll2"))),
  `loo2-with-rstan` = list(read_loo("real_wells.rds", c("log_lik_1", "log_lik_2"))),
  `loo2-weights` = milk_kline()
)

parallel_args <- function(mode, workers) {
  if (impl == "base") list(cores = if (mode == "serial") 1 else workers) else list(parallel = mode != "serial")
}

uses_waic <- c("bayesian-workflow/roaches", "loo2-weights")

vignette_code <- function(sets, mode, workers, with_waic) {
  pa <- parallel_args(mode, workers)
  for (models in sets) {
    loos <- lapply(models, function(m) {
      r_eff <- do.call(relative_eff, c(list(exp(m$ll)), if (!is.null(m$chain_id)) list(chain_id = m$chain_id), pa))
      suppressWarnings(do.call(loo, c(list(m$ll, r_eff = r_eff, save_psis = TRUE), pa)))
    })
    if (with_waic) lapply(models, function(m) suppressWarnings(do.call(waic, c(list(m$ll), pa))))
    if (length(loos) > 1) loo_compare(loos)
  }
  invisible(NULL)
}

repetitions <- 7
cat("impl,workload,mode,workers,median_s,p10_s,p90_s\n")
for (name in names(workloads)) {
  if (!is.na(args[4]) && name != args[4]) next
  for (setting in list(list("serial", 1), list("parallel", 2), list("parallel", 4), list("parallel", 8), list("parallel", 16))) {
    mode <- setting[[1]]
    workers <- setting[[2]]
    if (impl != "base" && mode == "parallel") {
      mirai::daemons(0)
      mirai::daemons(workers)
      invisible(mirai::call_mirai(mirai::everywhere(suppressMessages(library(loo)))))
    }
    with_waic <- name %in% uses_waic
    vignette_code(workloads[[name]], mode, workers, with_waic)
    elapsed <- replicate(repetitions, system.time(vignette_code(workloads[[name]], mode, workers, with_waic))[["elapsed"]])
    cat(sprintf("%s,%s,%s,%d,%.4f,%.4f,%.4f\n", impl, name, mode, workers, median(elapsed), quantile(elapsed, 0.1), quantile(elapsed, 0.9)))
    if (impl != "base" && mode == "parallel") mirai::daemons(0)
  }
}
