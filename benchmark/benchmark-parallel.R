# Benchmark loo's parallel code paths for a single installed version.
#
# Run once per version, pointing at an isolated library and tagging the output:
#   LOO_LIB=/tmp/loo-base-lib BENCH_LABEL=baseline Rscript benchmark/benchmark-parallel.R
#   LOO_LIB=/tmp/loo-new-lib  BENCH_LABEL=new      Rscript benchmark/benchmark-parallel.R
#
# Results are written to /tmp/bench-<label>.rds; compare.R aggregates them.
#
# The same user-facing calls (psis(), loo()) are timed for every version; the
# parallel backend (mclapply/parLapply vs mirai+mori) differs internally. For
# the new version we additionally time a "persist" mode where mirai::daemons()
# is started once before the timed iterations and reused (user-managed session
# pool), so the timed iterations measure steady-state cost with no per-call
# daemon spawn/teardown overhead.

lib <- Sys.getenv("LOO_LIB")
label <- Sys.getenv("BENCH_LABEL", unset = "unknown")
stopifnot(nzchar(lib))

suppressMessages({
  library(loo, lib.loc = lib)
  library(bench)
})

is_new <- identical(label, "new")
cores_grid <- c(1L, 4L, 6L)
iters <- 10L     # iterations for the small/cheap scenarios
big_iters <- 5L  # fewer iterations for the large, slow "worth parallelizing" scenarios

rows <- list()
record <- function(scenario, mode, cores, expr, n_iter = iters) {
  expr <- substitute(expr)
  pf <- parent.frame()
  eval(expr, pf) # warm up (process spawn / pool / JIT)
  b <- tryCatch(
    bench::mark(
      eval(expr, pf),
      iterations = n_iter, check = FALSE, memory = TRUE, filter_gc = FALSE
    ),
    error = function(e) {
      bench::mark(
        eval(expr, pf),
        iterations = n_iter, check = FALSE, memory = FALSE, filter_gc = FALSE
      )
    }
  )
  rows[[length(rows) + 1L]] <<- data.frame(
    label = label, scenario = scenario, mode = mode, cores = cores,
    iters = n_iter,
    median_s = as.numeric(b$median),
    mem_mb = if ("mem_alloc" %in% names(b) && !is.na(b$mem_alloc[1])) {
      as.numeric(b$mem_alloc) / 1e6
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  cat(sprintf(
    "  [%s] %-22s %-10s cores=%d  median=%.3fs\n",
    label, scenario, mode, cores, as.numeric(b$median)
  ))
}

# ---------------------------------------------------------------------------
# Scenario 1: standalone PSIS over a log-ratio matrix (partitioned columns).
#
# Counter-intuitively, the matrix `psis()` path is NOT worth parallelizing even
# when it is large: it has to ship a big log-ratio matrix out to the workers and
# return an equally large weighted matrix, so it is communication-bound. The
# parallel speedup stays ~1x no matter how big the problem is. The first two
# sizes are cheap exercisers; the last is large (~11s serial) and is included
# precisely to *demonstrate* that size alone does not make this path worth
# parallelizing -- compare its ~1x speedup with the loo.function scenarios below.
# ---------------------------------------------------------------------------
psis_sizes <- list(
  list(S = 2000, N = 1000, iters = iters),       # cheap exerciser
  list(S = 4000, N = 4000, iters = iters),       # cheap exerciser
  list(S = 6000, N = 10000, iters = big_iters)   # large but communication-bound (~1x)
)
for (sz in psis_sizes) {
  set.seed(2024) # identical inputs across versions
  S <- sz[["S"]]
  N <- sz[["N"]]
  n_iter <- sz[["iters"]]
  LL <- matrix(rnorm(S * N), nrow = S)
  re <- rep(1, N)
  scen <- sprintf("psis S=%d N=%d", S, N)
  for (k in cores_grid) {
    mode <- if (k == 1L) "serial" else "per-call"
    record(scen, mode, k, suppressWarnings(psis(-LL, r_eff = re, cores = k)), n_iter = n_iter)
    if (is_new && k > 1L) {
      mirai::daemons(k)
      record(scen, "persist", k, suppressWarnings(psis(-LL, r_eff = re, cores = k)), n_iter = n_iter)
      mirai::daemons(0)
    }
  }
  rm(LL)
  gc()
}

# ---------------------------------------------------------------------------
# Scenario 2: loo.function with a large broadcast `draws` matrix (the case
# where fork shares memory for free and mori must recover that benefit).
#
# This is the path that is genuinely worth parallelizing: `draws` is shared
# zero-copy across local workers via mori, and only tiny per-observation data
# and results move, so the per-observation work parallelizes cleanly. The first
# config is cheap (N=400, sub-second serial) and not worth parallelizing, but
# the two larger ones (~9s and ~18s serial) scale well even with the per-call
# pool -- ~2.3x at 4 cores and ~3x at 8 cores in local testing -- so a single
# one-off call already benefits.
# ---------------------------------------------------------------------------
llfun_b <- function(data_i, draws, ...) {
  dnorm(data_i$y, mean = draws[, "p1"], sd = abs(draws[, "p2"]) + 0.5, log = TRUE)
}
loo_sizes <- list(
  list(S = 8000, P = 150, N = 400, iters = iters),        # cheap: not worth parallelizing
  list(S = 8000, P = 200, N = 8000, iters = big_iters),   # ~9s serial: worth parallelizing
  list(S = 8000, P = 200, N = 10000, iters = big_iters)   # ~18s serial: worth parallelizing
)
for (sz in loo_sizes) {
  set.seed(7) # identical inputs across versions
  S2 <- sz[["S"]]
  P <- sz[["P"]]
  Nf <- sz[["N"]]
  n_iter <- sz[["iters"]]
  draws_big <- matrix(rnorm(S2 * P), nrow = S2, dimnames = list(NULL, paste0("p", seq_len(P))))
  data_f <- data.frame(y = rnorm(Nf))
  scen <- sprintf("loo.function S=%d P=%d N=%d (draws=%.1fMB)", S2, P, Nf, S2 * P * 8 / 1e6)
  for (k in cores_grid) {
    mode <- if (k == 1L) "serial" else "per-call"
    record(scen, mode, k, suppressWarnings(
      loo(llfun_b, data = data_f, draws = draws_big, cores = k)
    ), n_iter = n_iter)
    if (is_new && k > 1L) {
      mirai::daemons(k)
      record(scen, "persist", k, suppressWarnings(
        loo(llfun_b, data = data_f, draws = draws_big, cores = k)
      ), n_iter = n_iter)
      mirai::daemons(0)
    }
  }
  rm(draws_big)
  gc()
}

out <- do.call(rbind, rows)
saveRDS(out, sprintf("/tmp/bench-%s.rds", label))
cat(sprintf("\nSaved %d rows to /tmp/bench-%s.rds\n", nrow(out), label))
