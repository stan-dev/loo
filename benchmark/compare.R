# Aggregate and compare benchmark results from benchmark-parallel.R.
#   Rscript benchmark/compare.R
#
# Prints the comparison to the console and also writes a Markdown report
# (default /tmp/bench-comparison.md, override with the BENCH_MD env var) with
# the same two tables, for easy sharing.

# Directory containing this script, so the report lands next to it regardless
# of the working directory the script is launched from.
.this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(.this_file) == 1L) {
  dirname(normalizePath(.this_file))
} else {
  "benchmark"
}

md_out <- Sys.getenv("BENCH_MD", unset = file.path(script_dir, "bench-comparison.md"))
peak_out <- Sys.getenv("PEAK_OUT", unset = "/tmp/bench-peakmem.tsv")

base <- readRDS("/tmp/bench-baseline.rds")
new <- readRDS("/tmp/bench-new.rds")
all <- rbind(base, new)

key <- function(d) paste(d$scenario, d$cores)

# Median seconds, keyed by version/mode.
med <- function(df, lab, mode) {
  sel <- df[df$label == lab & df$mode == mode, ]
  setNames(sel$median_s, key(sel))
}
mem <- function(df, lab, mode) {
  sel <- df[df$label == lab & df$mode == mode, ]
  setNames(sel$mem_mb, key(sel))
}

scen_cores <- unique(all[, c("scenario", "cores")])
scen_cores <- scen_cores[order(scen_cores$scenario, scen_cores$cores), ]

# Iterations per measurement (recorded by benchmark-parallel.R). Older result
# files may predate this column.
iters_used <- if ("iters" %in% names(all)) {
  iv <- sort(unique(all$iters[!is.na(all$iters)]))
  if (length(iv) == 1L) as.character(iv) else paste(range(iv), collapse = "-")
} else {
  NA_character_
}
iters_label <- if (is.na(iters_used)) "a few" else iters_used

b_serial <- med(base, "baseline", "serial")
b_par <- med(base, "baseline", "per-call")
n_serial <- med(new, "new", "serial")
n_call <- med(new, "new", "per-call")
n_persist <- med(new, "new", "persist")

base_med <- c(b_serial, b_par)
new_call_med <- c(n_serial, n_call)

get1 <- function(vec, k) {
  if (!is.null(vec) && k %in% names(vec)) vec[[k]] else NA_real_
}
fmt <- function(x) ifelse(is.na(x), "      -", sprintf("%7.3f", x))
spd <- function(num, den) ifelse(is.na(num) | is.na(den), "    -", sprintf("%4.2fx", num / den))

cat(sprintf("\n=== Median wall-clock time (s) and speedup vs baseline (%s iterations) ===\n", iters_label))
cat(sprintf(
  "%-42s %5s | %8s %8s %8s | %7s %7s\n",
  "scenario", "cores", "base", "new/call", "new/per",
  "call", "persist"
))
cat(strrep("-", 96), "\n")
for (i in seq_len(nrow(scen_cores))) {
  k <- paste(scen_cores$scenario[i], scen_cores$cores[i])
  bm <- get1(base_med, k)
  nc <- get1(new_call_med, k)
  np <- get1(n_persist, k)
  cat(sprintf(
    "%-42s %5d | %8s %8s %8s | %7s %7s\n",
    scen_cores$scenario[i], scen_cores$cores[i],
    fmt(bm), fmt(nc), fmt(np),
    spd(bm, nc), spd(bm, np)
  ))
}

cat("\n=== Main-process memory allocation (MB) ===\n")
b_mem <- c(mem(base, "baseline", "serial"), mem(base, "baseline", "per-call"))
n_mem <- c(mem(new, "new", "serial"), mem(new, "new", "per-call"))
p_mem <- mem(new, "new", "persist")
cat(sprintf("%-42s %5s | %8s %8s %8s\n", "scenario", "cores", "base", "new/call", "new/per"))
cat(strrep("-", 80), "\n")
for (i in seq_len(nrow(scen_cores))) {
  k <- paste(scen_cores$scenario[i], scen_cores$cores[i])
  cat(sprintf(
    "%-42s %5d | %8s %8s %8s\n",
    scen_cores$scenario[i], scen_cores$cores[i],
    fmt(get1(b_mem, k)), fmt(get1(n_mem, k)), fmt(get1(p_mem, k))
  ))
}
cat("\nspeedup > 1 means new is faster than baseline.\n")


# Peak-memory results (optional) ---------------------------------------------
# peak-mem.sh appends "<label>\t<mode>\t<cores>\t<peak_mb>" rows to peak_out.

peak <- NULL
if (file.exists(peak_out) && file.info(peak_out)$size > 0) {
  peak <- utils::read.delim(
    peak_out, header = FALSE, stringsAsFactors = FALSE,
    col.names = c("label", "mode", "cores", "peak_mb")
  )
  # Keep only the most recent record for each label/mode/cores combination.
  peak$.k <- paste(peak$label, peak$mode, peak$cores)
  peak <- peak[!duplicated(peak$.k, fromLast = TRUE), ]
  peak <- peak[order(peak$cores, peak$label, peak$mode), ]

  cat("\n=== Peak RSS of whole process tree (MB) ===\n")
  cat(sprintf("%-10s %-10s %5s | %10s\n", "label", "mode", "cores", "peak_mb"))
  cat(strrep("-", 44), "\n")
  for (i in seq_len(nrow(peak))) {
    cat(sprintf(
      "%-10s %-10s %5s | %10.0f\n",
      peak$label[i], peak$mode[i], peak$cores[i], peak$peak_mb[i]
    ))
  }
} else {
  cat(sprintf(
    "\n(no peak-memory results found at %s; run benchmark/peak-mem.sh to add them)\n",
    peak_out
  ))
}


# Markdown report ------------------------------------------------------------
# Same numbers as above, formatted as Markdown tables for easy sharing.

md_num <- function(x) ifelse(is.na(x), "—", sprintf("%.3f", x))
md_int <- function(x) ifelse(is.na(x), "—", sprintf("%.0f", x))
md_spd <- function(num, den) {
  ifelse(is.na(num) | is.na(den), "—", sprintf("%.2fx", num / den))
}

md <- c(
  "# loo parallel benchmark comparison",
  "",
  sprintf("_Generated %s._", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "**Columns / modes.** `base` is the baseline version. `new/call` is the new version's default per-call pool (created and torn down each call). `new/persist` is the new version with a user-managed `mirai::daemons()` session pool, reused across calls. `cores = 1` rows are fully serial (the parallel backend is never used).",
  "",
  "## Median wall-clock time (s) and speedup vs baseline",
  "",
  sprintf("Median over %s iterations of one `psis()`/`loo()` call (a warm-up run is excluded). `speedup = base / new`, so a value `> 1` means the new version is faster. Expect `new/persist` to win when many calls reuse the pool, and `new/call` to look slower than `base` for cheap problems because it pays pool start-up/teardown on every call. The `cores = 1` `base` and `new/call` numbers should be roughly equal (both serial); sizeable gaps there are run-to-run noise, not real differences.", iters_label),
  "",
  "| scenario | cores | base | new/call | new/persist | speedup (call) | speedup (persist) |",
  "|:---|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(scen_cores))) {
  k <- paste(scen_cores$scenario[i], scen_cores$cores[i])
  bm <- get1(base_med, k)
  nc <- get1(new_call_med, k)
  np <- get1(n_persist, k)
  md <- c(md, sprintf(
    "| %s | %d | %s | %s | %s | %s | %s |",
    scen_cores$scenario[i], scen_cores$cores[i],
    md_num(bm), md_num(nc), md_num(np),
    md_spd(bm, nc), md_spd(bm, np)
  ))
}

md <- c(
  md,
  "",
  "## Main-process memory allocation (MB)",
  "",
  "Total bytes allocated on the R heap by the *coordinator* process during the call (cumulative churn, **not** peak and **not** net), as measured by `bench`'s allocation profiler. It does **not** include memory used by worker processes, nor off-heap memory such as the `mori` shared-memory segment for the broadcast `draws`. That is why parallel rows are tiny: the heavy allocation happens in the workers, out of the profiler's view. Use it to gauge allocation pressure on the main process, not total footprint.",
  "",
  "| scenario | cores | base | new/call | new/persist |",
  "|:---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(scen_cores))) {
  k <- paste(scen_cores$scenario[i], scen_cores$cores[i])
  md <- c(md, sprintf(
    "| %s | %d | %s | %s | %s |",
    scen_cores$scenario[i], scen_cores$cores[i],
    md_num(get1(b_mem, k)), md_num(get1(n_mem, k)), md_num(get1(p_mem, k))
  ))
}

# Peak-RSS table (only when peak-mem.sh results are available).
if (!is.null(peak)) {
  md <- c(
    md,
    "",
    "## Peak RSS of the whole process tree (MB)",
    "",
    "From `peak-mem.sh` (single large-`draws` `loo()` run; Linux only). Maximum summed resident memory of the *entire* process tree (main process plus all workers), sampled during the run. This is the metric for the job's real memory footprint. On a local pool, `mori` shares the `draws` matrix across workers (zero-copy), so peak RSS stays close to a single copy rather than growing with the number of workers.",
    "",
    "| label | mode | cores | peak RSS (MB) |",
    "|:---|:---|---:|---:|"
  )
  for (i in seq_len(nrow(peak))) {
    md <- c(md, sprintf(
      "| %s | %s | %s | %s |",
      peak$label[i], peak$mode[i], peak$cores[i], md_int(peak$peak_mb[i])
    ))
  }
}

# Caveats footer.
md <- c(
  md,
  "",
  "## Caveats",
  "",
  "Results are platform-dependent and sensitive to machine load; run the baseline and new versions back to back on an idle machine, and ignore differences smaller than the run-to-run noise.",
  ""
)

writeLines(md, md_out)
cat(sprintf("\nMarkdown report written to %s\n", md_out))
