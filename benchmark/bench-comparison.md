# loo parallel benchmark comparison

_Generated 2026-07-01 11:23:20._

**Columns / modes.** `base` is the baseline version. `new/call` is the new version's default per-call pool (created and torn down each call). `new/persist` is the new version with a user-managed `mirai::daemons()` session pool, reused across calls. `cores = 1` rows are fully serial (the parallel backend is never used).

## Median wall-clock time (s) and speedup vs baseline

Median over 5-10 iterations of one `psis()`/`loo()` call (a warm-up run is excluded). `speedup = base / new`, so a value `> 1` means the new version is faster. Expect `new/persist` to win when many calls reuse the pool, and `new/call` to look slower than `base` for cheap problems because it pays pool start-up/teardown on every call. The `cores = 1` `base` and `new/call` numbers should be roughly equal (both serial); sizeable gaps there are run-to-run noise, not real differences.

| scenario | cores | base | new/call | new/persist | speedup (call) | speedup (persist) |
|:---|---:|---:|---:|---:|---:|---:|
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 1 | 0.677 | 0.678 | — | 1.00x | — |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 4 | 0.468 | 1.161 | 0.231 | 0.40x | 2.02x |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 6 | 0.379 | 1.200 | 0.147 | 0.32x | 2.58x |
| loo.function S=8000 P=200 N=10000 (draws=12.8MB) | 1 | 15.150 | 17.807 | — | 0.85x | — |
| loo.function S=8000 P=200 N=10000 (draws=12.8MB) | 4 | 4.716 | 5.691 | 5.013 | 0.83x | 0.94x |
| loo.function S=8000 P=200 N=10000 (draws=12.8MB) | 6 | 2.598 | 4.036 | 3.258 | 0.64x | 0.80x |
| loo.function S=8000 P=200 N=8000 (draws=12.8MB) | 1 | 13.527 | 13.847 | — | 0.98x | — |
| loo.function S=8000 P=200 N=8000 (draws=12.8MB) | 4 | 4.325 | 4.944 | 4.086 | 0.87x | 1.06x |
| loo.function S=8000 P=200 N=8000 (draws=12.8MB) | 6 | 2.597 | 3.437 | 2.650 | 0.76x | 0.98x |
| psis S=2000 N=1000 | 1 | 0.388 | 0.371 | — | 1.05x | — |
| psis S=2000 N=1000 | 4 | 0.242 | 1.153 | 0.214 | 0.21x | 1.13x |
| psis S=2000 N=1000 | 6 | 0.252 | 1.190 | 0.176 | 0.21x | 1.43x |
| psis S=4000 N=4000 | 1 | 2.892 | 2.764 | — | 1.05x | — |
| psis S=4000 N=4000 | 4 | 1.704 | 2.331 | 1.465 | 0.73x | 1.16x |
| psis S=4000 N=4000 | 6 | 1.367 | 2.547 | 1.786 | 0.54x | 0.77x |
| psis S=6000 N=10000 | 1 | 9.641 | 9.195 | — | 1.05x | — |
| psis S=6000 N=10000 | 4 | 6.203 | 7.735 | 7.242 | 0.80x | 0.86x |
| psis S=6000 N=10000 | 6 | 5.384 | 7.408 | 6.801 | 0.73x | 0.79x |

## Main-process memory allocation (MB)

Total bytes allocated on the R heap by the *coordinator* process during the call (cumulative churn, **not** peak and **not** net), as measured by `bench`'s allocation profiler. It does **not** include memory used by worker processes, nor off-heap memory such as the `mori` shared-memory segment for the broadcast `draws`. That is why parallel rows are tiny: the heavy allocation happens in the workers, out of the profiler's view. Use it to gauge allocation pressure on the main process, not total footprint.

| scenario | cores | base | new/call | new/persist |
|:---|---:|---:|---:|---:|
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 1 | 803.292 | 803.292 | — |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 4 | — | 1.112 | 1.092 |
| loo.function S=8000 P=150 N=400 (draws=9.6MB) | 6 | — | 1.596 | 1.576 |
| loo.function S=8000 P=200 N=10000 (draws=12.8MB) | 1 | 20082.281 | 20082.281 | — |
| loo.function S=8000 P=200 N=10000 (draws=12.8MB) | 4 | — | 2.917 | 2.897 |
| loo.function S=8000 P=200 N=10000 (draws=12.8MB) | 6 | — | 3.401 | 3.381 |
| loo.function S=8000 P=200 N=8000 (draws=12.8MB) | 1 | 16065.825 | 16065.825 | — |
| loo.function S=8000 P=200 N=8000 (draws=12.8MB) | 4 | — | 2.541 | 2.521 |
| loo.function S=8000 P=200 N=8000 (draws=12.8MB) | 6 | — | 3.025 | 3.005 |
| psis S=2000 N=1000 | 1 | 304.130 | 304.130 | — |
| psis S=2000 N=1000 | 4 | — | 120.986 | 120.966 |
| psis S=2000 N=1000 | 6 | — | 121.335 | 121.315 |
| psis S=4000 N=4000 | 1 | 2306.371 | 2306.371 | — |
| psis S=4000 N=4000 | 4 | — | 961.700 | 961.679 |
| psis S=4000 N=4000 | 6 | — | 962.048 | 962.028 |
| psis S=6000 N=10000 | 1 | 8458.151 | 8458.151 | — |
| psis S=6000 N=10000 | 4 | — | 3603.208 | 3603.188 |
| psis S=6000 N=10000 | 6 | — | 3603.557 | 3603.537 |

## Caveats

Results are platform-dependent and sensitive to machine load; run the baseline and new versions back to back on an idle machine, and ignore differences smaller than the run-to-run noise.

