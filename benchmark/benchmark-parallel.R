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
# the new version we additionally time a "persist" mode that opts in to loo's
# persistent session pool via `options(loo.daemons = k)` (equivalently the
# `LOO_DAEMONS` environment variable). loo then creates the local mirai pool
# lazily on the first (warm-up) call and reuses it for every later call, so the
# timed iterations measure steady-state cost with no per-call daemon
# spawn/teardown overhead.

lib <- Sys.getenv("LOO_LIB")
label <- Sys.getenv("BENCH_LABEL", unset = "unknown")
stopifnot(nzchar(lib))

suppressMessages({
  library(loo, lib.loc = lib)
  library(bench)
})

is_new <- identical(label, "new")
cores_grid <- c(1L, 4L, 8L)
iters <- 10L

rows <- list()
record <- function(scenario, mode, cores, expr) {
  expr <- substitute(expr)
  pf <- parent.frame()
  eval(expr, pf) # warm up (process spawn / pool / JIT)
  b <- tryCatch(
    bench::mark(
      eval(expr, pf),
      iterations = iters, check = FALSE, memory = TRUE, filter_gc = FALSE
    ),
    error = function(e) {
      bench::mark(
        eval(expr, pf),
        iterations = iters, check = FALSE, memory = FALSE, filter_gc = FALSE
      )
    }
  )
  rows[[length(rows) + 1L]] <<- data.frame(
    label = label, scenario = scenario, mode = mode, cores = cores,
    iters = iters,
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
# ---------------------------------------------------------------------------
psis_sizes <- list(c(S = 2000, N = 1000), c(S = 4000, N = 4000))
for (sz in psis_sizes) {
  set.seed(2024) # identical inputs across versions
  S <- sz[["S"]]
  N <- sz[["N"]]
  LL <- matrix(rnorm(S * N), nrow = S)
  re <- rep(1, N)
  scen <- sprintf("psis S=%d N=%d", S, N)
  for (k in cores_grid) {
    mode <- if (k == 1L) "serial" else "per-call"
    record(scen, mode, k, suppressWarnings(psis(-LL, r_eff = re, cores = k)))
    if (is_new && k > 1L) {
      # Opt in to loo's persistent session pool; the warm-up call inside
      # record() creates it and the timed iterations reuse it.
      old_opt <- options(loo.daemons = k)
      record(scen, "persist", k, suppressWarnings(psis(-LL, r_eff = re, cores = k)))
      mirai::daemons(0)
      options(old_opt)
    }
  }
}

# ---------------------------------------------------------------------------
# Scenario 2: loo.function with a large broadcast `draws` matrix (the case
# where fork shares memory for free and mori must recover that benefit).
# ---------------------------------------------------------------------------
set.seed(7)
S2 <- 8000L
P <- 150L
Nf <- 400L
draws_big <- matrix(rnorm(S2 * P), nrow = S2, dimnames = list(NULL, paste0("p", seq_len(P))))
data_f <- data.frame(y = rnorm(Nf))
llfun_b <- function(data_i, draws, ...) {
  dnorm(data_i$y, mean = draws[, "p1"], sd = abs(draws[, "p2"]) + 0.5, log = TRUE)
}
scen <- sprintf("loo.function S=%d P=%d N=%d (draws=%.1fMB)", S2, P, Nf, S2 * P * 8 / 1e6)
for (k in cores_grid) {
  mode <- if (k == 1L) "serial" else "per-call"
  record(scen, mode, k, suppressWarnings(
    loo(llfun_b, data = data_f, draws = draws_big, cores = k)
  ))
  if (is_new && k > 1L) {
    # Opt in to loo's persistent session pool; the warm-up call inside
    # record() creates it and the timed iterations reuse it.
    old_opt <- options(loo.daemons = k)
    record(scen, "persist", k, suppressWarnings(
      loo(llfun_b, data = data_f, draws = draws_big, cores = k)
    ))
    mirai::daemons(0)
    options(old_opt)
  }
}

out <- do.call(rbind, rows)
saveRDS(out, sprintf("/tmp/bench-%s.rds", label))
cat(sprintf("\nSaved %d rows to /tmp/bench-%s.rds\n", nrow(out), label))
