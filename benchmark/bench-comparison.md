# loo parallel benchmark comparison

_Generated 2026-06-30 16:29:48._

**Columns / modes.** `base` is the baseline version. `new/call` is the new version's default per-call pool (created and torn down each call). `new/persist` is the new version's persistent session pool (`options(loo.daemons = k)`), reused across calls. `cores = 1` rows are fully serial (the parallel backend is never used).

## Median wall-clock time (s) and speedup vs baseline

Median over 10 iterations of one `psis()`/`loo()` call (a warm-up run is excluded). `speedup = base / new`, so a value `> 1` means the new version is faster. Expect `new/persist` to win when many calls reuse the pool, and `new/call` to look slower than `base` for cheap problems because it pays pool start-up/teardown on every call. The `cores = 1` `base` and `new/call` numbers should be roughly equal (both serial); sizeable gaps there are run-to-run noise, not real differences.

| scenario | cores | base | new/call | new/persist | speedup (call) | speedup (persist) |
|:---|---:|---:|---:|---:|---:|---:|
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 1 | 0.458 | 0.406 | — | 1.13x | — |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 4 | 0.670 | 1.315 | 0.175 | 0.51x | 3.83x |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 8 | 0.536 | 1.433 | 0.133 | 0.37x | 4.02x |
| psis S=2000 N=1000 | 1 | 0.416 | 0.403 | — | 1.03x | — |
| psis S=2000 N=1000 | 4 | 0.237 | 1.262 | 0.158 | 0.19x | 1.51x |
| psis S=2000 N=1000 | 8 | 0.205 | 1.276 | 0.136 | 0.16x | 1.51x |
| psis S=4000 N=4000 | 1 | 2.998 | 2.628 | — | 1.14x | — |
| psis S=4000 N=4000 | 4 | 2.253 | 2.385 | 1.455 | 0.94x | 1.55x |
| psis S=4000 N=4000 | 8 | 1.707 | 2.258 | 1.050 | 0.76x | 1.63x |

## Main-process memory allocation (MB)

Total bytes allocated on the R heap by the *coordinator* process during the call (cumulative churn, **not** peak and **not** net), as measured by `bench`'s allocation profiler. It does **not** include memory used by worker processes, nor off-heap memory such as the `mori` shared-memory segment for the broadcast `draws`. That is why parallel rows are tiny: the heavy allocation happens in the workers, out of the profiler's view. Use it to gauge allocation pressure on the main process, not total footprint.

| scenario | cores | base | new/call | new/persist |
|:---|---:|---:|---:|---:|
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 1 | 803.292 | 803.292 | — |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 4 | — | 0.720 | 0.700 |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 8 | — | 1.297 | 1.276 |
| psis S=2000 N=1000 | 1 | 304.130 | 304.130 | — |
| psis S=2000 N=1000 | 4 | — | 120.805 | 120.785 |
| psis S=2000 N=1000 | 8 | — | 121.322 | 121.302 |
| psis S=4000 N=4000 | 1 | 2306.371 | 2306.371 | — |
| psis S=4000 N=4000 | 4 | — | 961.519 | 961.498 |
| psis S=4000 N=4000 | 8 | — | 962.035 | 962.015 |

## Peak RSS of the whole process tree (MB)

From `peak-mem.sh` (single large-`draws` `loo()` run; Linux only). Maximum summed resident memory of the *entire* process tree (main process plus all workers), sampled during the run. This is the metric for the job's real memory footprint. On a local pool, `mori` shares the `draws` matrix across workers (zero-copy), so peak RSS stays close to a single copy rather than growing with the number of workers.

| label | mode | cores | peak RSS (MB) |
|:---|:---|---:|---:|
| baseline | per-call | 8 | 2418 |
| new | per-call | 8 | 343 |
| new | persist | 8 | 343 |

## Caveats

Results are platform-dependent and sensitive to machine load; run the baseline and new versions back to back on an idle machine, and ignore differences smaller than the run-to-run noise.

