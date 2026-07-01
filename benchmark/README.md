# loo parallel benchmarks

These scripts measure the performance of loo's parallel code paths and compare
two installed versions of the package side by side:

- **`baseline`** — a pre-`mirai` version (the old `mclapply`/`parLapply`
  backend), e.g. the released version from CRAN.
- **`new`** — the current working tree (the `mirai` + `mori` backend, including
  the persistent session pool controlled by `options(loo.daemons = k)` /
  `LOO_DAEMONS`).

The same user-facing calls (`psis()`, `loo()`) are timed for every version; only
the internal parallel backend differs. For the `new` version we additionally
time a **persist** mode that opts in to the persistent session pool, so we can
separate per-call daemon spawn/teardown overhead from the steady-state cost.

## Files

| File | Purpose |
|---|---|
| `benchmark-parallel.R` | Times `psis()`/`loo()` across cores for one installed version; writes `/tmp/bench-<label>.rds`. |
| `compare.R` | Reads the two `.rds` files and prints a wall-clock + memory comparison table. |
| `peak-mem-run.R` | A single large-`draws` `loo()` run used for peak-memory measurement. |
| `peak-mem.sh` | Linux-only driver that samples the RSS of the whole process tree during `peak-mem-run.R`. |

## Prerequisites

The benchmark scripts need the `bench` package (in addition to whatever loo
needs):

```bash
Rscript -e 'install.packages("bench")'
```

## Step 1 — Install each version into its own library

Each run reads the library path from the `LOO_LIB` environment variable, so put
each version in its own directory:

```bash
mkdir -p /tmp/loo-base-lib /tmp/loo-new-lib

# "new" = this working tree
R CMD INSTALL --library=/tmp/loo-new-lib .

# "baseline" = a pre-mirai version to compare against (CRAN release shown here;
# alternatively check out an older git ref elsewhere and install that)
Rscript -e 'install.packages("loo", lib = "/tmp/loo-base-lib")'
```

## Step 2 — Run the timing benchmark once per version

Run from the package root. `BENCH_LABEL` tags the output file
(`/tmp/bench-<label>.rds`):

```bash
LOO_LIB=/tmp/loo-base-lib BENCH_LABEL=baseline Rscript benchmark/benchmark-parallel.R
LOO_LIB=/tmp/loo-new-lib  BENCH_LABEL=new      Rscript benchmark/benchmark-parallel.R
```

The `persist` mode is only measured for `BENCH_LABEL=new`, because
`options(loo.daemons)` is a no-op in the baseline version.

## Step 3 — Aggregate and compare

```bash
Rscript benchmark/compare.R
```

This prints two tables (median wall-clock time with speedups, and main-process
memory) with these columns:

- `base` — baseline version.
- `new/call` — new version, default per-call pool (created and torn down each
  call).
- `new/per` — new version, persistent session pool (`options(loo.daemons = k)`),
  reused across calls.

A speedup `> 1` means the new version is faster than the baseline. The report
also states how many iterations each median is based on (recorded from
`benchmark-parallel.R`'s `iters` setting).

In addition to the console output, `compare.R` writes a Markdown report with the
same tables (handy for pasting into issues/PRs). By default it is written next to
the script as `benchmark/bench-comparison.md`; override the path with the
`BENCH_MD` environment variable:

```bash
BENCH_MD=/tmp/benchmark-results.md Rscript benchmark/compare.R
```

## Optional — Peak memory (Linux only)

`peak-mem.sh` samples the total RSS of the R process and all of its workers
during one large-`draws` `loo()` run. It takes positional arguments
`LOO_LIB BENCH_LABEL MODE CORES`, where `MODE` is `per-call` or `persist`:

```bash
benchmark/peak-mem.sh /tmp/loo-base-lib baseline per-call 8
benchmark/peak-mem.sh /tmp/loo-new-lib  new      per-call 8
benchmark/peak-mem.sh /tmp/loo-new-lib  new      persist  8
```

For the `new` version it also prints the `mori` transport size of the broadcast
`draws` object (raw serialized MB vs the shared-memory reference in bytes),
which shows the zero-copy benefit on a local pool.

Each run also appends its peak-RSS result to a tab-separated file (default
`/tmp/bench-peakmem.tsv`, override with `PEAK_OUT`). If that file is present when
you run `compare.R`, the peak-RSS numbers are folded into the report as an extra
"Peak RSS of the whole process tree" table. Delete the file between fresh runs so
stale rows aren't mixed in (only the most recent row per label/mode/cores is
kept).

The Markdown report also ends with a **"How to read these numbers"** section
explaining the difference between median time/speedup, main-process heap
allocation (churn, not peak, workers excluded), and peak RSS (the whole process
tree's real footprint).

## Tuning

Edit the top of `benchmark-parallel.R` to match your machine / problem sizes:

- `cores_grid` (default `c(1, 4, 6)`) — the core counts to sweep.
- `iters` (default `10`) — iterations per `bench::mark()` measurement for the
  small/cheap scenarios.
- `big_iters` (default `5`) — iterations for the large, slow scenarios that are
  actually worth parallelizing (fewer iterations keeps total runtime sane).
- `psis_sizes` and `loo_sizes` — the problem sizes. Each entry carries its own
  `iters`.

Which scenarios are actually worth parallelizing:

- **`loo.function` (`loo_sizes`)** is the path that parallelizes well. `draws`
  is shared zero-copy across local workers via `mori`, so only tiny
  per-observation data/results move. The two larger configs (~9s and ~18s
  serial) reach roughly 2.3x at 4 cores and 3x at 8 cores in local testing —
  even with the per-call pool, so a single one-off call already benefits.
- **Matrix `psis` (`psis_sizes`)** is *not* worth parallelizing, even when
  large: it must ship a big log-ratio matrix out to the workers and return an
  equally large weighted matrix, so it is communication-bound and stays near 1x
  regardless of size. The large `S=6000 N=20000` entry (~960 MB) is included
  deliberately to show this — make sure the machine has enough RAM before adding
  bigger ones.

For `peak-mem-run.R`, adjust `S`, `P`, and `Nf` to change the size of the
broadcast `draws` matrix.

## Notes

- All output is written to `/tmp`. Remove `/tmp/bench-*.rds` between experiments
  if you change the scenarios, so stale rows aren't mixed in.
- Results are platform-dependent; run baseline and new back to back on an idle
  machine for a fair comparison.
