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
