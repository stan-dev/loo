# Single large-draws loo.function run for peak-memory measurement.
# Driven by peak-mem.sh, which samples the total RSS of this process tree.
#   LOO_LIB=... BENCH_LABEL=baseline|new MODE=per-call|persist CORES=8 Rscript peak-mem-run.R
#
# In "persist" mode (new version only) mirai::daemons() is started before the
# loo() call so the pool is reused for the duration of the run.

lib <- Sys.getenv("LOO_LIB")
label <- Sys.getenv("BENCH_LABEL", unset = "unknown")
mode <- Sys.getenv("MODE", unset = "per-call")
cores <- as.integer(Sys.getenv("CORES", unset = "8"))
stopifnot(nzchar(lib))
suppressMessages(library(loo, lib.loc = lib))

set.seed(7)
S <- 30000L
P <- 400L # draws ~ 30000 * 400 * 8 = 96 MB
Nf <- 300L
draws_big <- matrix(rnorm(S * P), nrow = S, dimnames = list(NULL, paste0("p", seq_len(P))))
data_f <- data.frame(y = rnorm(Nf))
llfun_b <- function(data_i, draws, ...) {
  dnorm(data_i$y, mean = draws[, "p1"], sd = abs(draws[, "p2"]) + 0.5, log = TRUE)
}

# Per-worker transport size of the broadcast object (mori metric).
if (label == "new") {
  raw_bytes <- length(serialize(draws_big, NULL))
  shared <- mori::share(draws_big)
  ref_bytes <- length(serialize(shared, NULL))
  cat(sprintf(
    "TRANSPORT draws raw=%.1fMB  mori_ref=%d bytes  (%.0fx smaller)\n",
    raw_bytes / 1e6, ref_bytes, raw_bytes / ref_bytes
  ))
}

if (label == "new" && mode == "persist") {
  mirai::daemons(cores)
}

invisible(suppressWarnings(loo(llfun_b, data = data_f, draws = draws_big, cores = cores)))

if (label == "new" && mode == "persist") {
  mirai::daemons(0)
}
cat("RUN COMPLETE\n")
