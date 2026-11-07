# Benchmark checks

Experiments run while building `bench/run.sh`. Windows 11, i7-12700H, plugged in, 4 workers. Not part of the PR results.

## Function-based WAIC: old `cores` code against the new `loo_map`

`waic(llfun, data, draws)`, 4000 observations, 2000 draws, fresh R process per row, three rounds. Milliseconds.

| | Start-up | Warm median per call |
|---|---|---|
| v3 `cores = 4` (PSOCK cluster built on every call) | none | 760, 804, 1583 |
| new, mirai | 390 to 730 | 491, 528, 908 |
| new, mizu | 394 to 850 | 521, 541, 885 |
| sequential | none | 1108, 1408, 2054 |

The same sequential call varied by almost 2x between runs, so single timings on this laptop are noisy.

## Batching variants that keep warnings

Median ms over 11 interleaved rounds (variant order shuffled each round). Every variant matched the sequential result and re-signalled all warnings.

| Workload | mirai blocks (current) | mirai chunks of 250 | mirai chunks of 50 | mizu per element | mizu blocks | mirai per element |
|---|---|---|---|---|---|---|
| WAIC 4000 x 2000 | 969 | 972 | 915 | 1193 | 1138 | 2427 |
| WAIC 20000 x 200 | 963 | 956 | 1186 | 1209 | 1003 | not run |
| `loo()` function 2000 x 1000 | 2095 | 1950 | 2097 | 2133 | 3622 | 2859 |

The IQRs overlap for the top rows, so none beats the current one block per daemon.

## Repetitions

One 40-repetition run per benchmark, re-analyzed on its first n repetitions. Width of the 80% interval of the relative change, in percentage points. Later repetitions ran while other jobs used the machine.

| Benchmark | Metric | n=10 | n=20 | n=30 | n=40 |
|---|---|---|---|---|---|
| `vig_milk_kline` | time | 0.36 | 0.82 | 0.71 | 0.59 |
| | memory | 0.44 | 0.57 | 0.51 | 0.43 |
| `relative_eff` | time | 7.13 | 6.14 | 5.84 | 6.16 |
| | memory | 0.49 | 0.45 | 0.32 | 0.28 |
| `loo_function` | time | 10.0 | 13.5 | 14.6 | 14.6 |
| | memory | 0.33 | 0.24 | 0.21 | 0.19 |

Memory estimates are tight and stable from 10 repetitions. Time intervals are set by machine noise more than by n.

## Sleep ruins runs

In a mori test one repetition took 8,867 s (candidate) and 55 s (baseline) while Windows logged connected-standby events. `bench/run.sh` now holds the machine awake while it runs (`SetThreadExecutionState` on Windows, `caffeinate` on macOS, `systemd-inhibit` on Linux).

## Data files

The `real_*.rds` files in `bench/` total 551 MB. `qs2` with zstd level 19 and byte shuffle gives 476 MB (13% smaller); wells and wells4 stay at 153 MB each, over GitHub's 100 MB file limit. They are MCMC output and do not compress losslessly.

## Linux slowdowns are daemon start-up, not per-call cost
foil times the whole benchmark script, and each candidate script starts its daemons inside it, while the baseline forks for free on Linux. On the cluster the candidate takes 1.2 to 1.7 s on the short benchmarks whatever the worker count, which is the added start-up.

Checked on WSL (Ubuntu, R 4.5.3, 20 threads) with `waic()` on the function interface (4000 observations, 2000 draws), median of 7 calls, daemons already running:

| Workers | Baseline `cores = N` (fork) | Candidate `parallel = TRUE` |
|---|---|---|
| 1 | 0.324 s | 0.337 s (`parallel = FALSE` path) |
| 2 | 0.209 s | 0.170 s |
| 8 | 0.146 s | 0.069 s |
| 16 | 0.123 s | 0.054 s |

Starting daemons and loading loo on them took 0.33 s for 2, 0.37 s for 8 and 0.74 s for 24 on the same machine. Per call the candidate is faster than forking; the benchmark scripts pay the start-up once or twice per run.

Closure size and the other function-based paths, same WSL setup, `loo()` and `loo_subsample()` (300 observations, `plpd`) on the function interface, 4000 observations, median of 5 calls with daemons running:

| Workers | `loo()` baseline | `loo()` candidate | `loo_subsample()` baseline | `loo_subsample()` candidate |
|---|---|---|---|---|
| 2 | 1.750 s | 1.471 s | 0.548 s | 0.172 s |
| 8 | 1.187 s | 0.719 s | 0.585 s | 0.077 s |

The log-likelihood function serializes to 412 bytes (data 64 KB, draws 32 KB), so nothing large rides along with the closure. Benchmark scripts define it in the global environment; a function created inside another function would carry that frame's objects.
