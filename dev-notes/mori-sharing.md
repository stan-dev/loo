# When loo_map shares arguments with mori

`loo_map` shares an argument with `mori::share()` when daemons are local and either the tasks are sent one per element (`batch = FALSE`) or the argument is larger than `share_bytes` bytes. Smaller arguments are copied. `share_bytes` defaults to `getOption("loo.share_bytes", 4 * 2^20)`, which is 4 MiB. `options(loo.parallel = ...)` sets the default of `parallel` in the exported functions.

## Measurements

750 `llfun`-style subsetting calls in one daemon, and sending the object to 4 daemons. Seconds.

| Matrix | Windows compute raw / shared | Windows send raw / shared | WSL compute raw / shared | WSL send raw / shared |
|---|---|---|---|---|
| 0.03 MiB | 0.030 / 0.050 | 0.000 / 0.000 | 0.013 / 0.047 | 0.002 / 0.002 |
| 0.31 MiB | 0.030 / 0.060 | 0.000 / 0.000 | 0.013 / 0.052 | 0.005 / 0.002 |
| 3 MiB | 0.050 / 0.110 | 0.020 / 0.000 | 0.029 / 0.105 | 0.035 / 0.002 |
| 15 MiB | 0.080 / 0.130 | 0.140 / 0.000 | 0.054 / 0.098 | 0.158 / 0.002 |
| 61 MiB | 0.150 / 0.140 | 0.480 / 0.000 | 0.148 / 0.120 | 0.558 / 0.002 |
| 153 MiB | 0.300 / 0.200 | 1.550 / 0.000 | 0.247 / 0.142 | 1.357 / 0.002 |

Reading a shared object costs a roughly fixed 40 to 75 ms per 750 accesses, at any size. Copying costs grow with size. They cross between 3 and 15 MiB.

One task per element, 4 daemons, 5 ms of work per task (Windows): copy and share take the same time for arguments up to 2 MiB at 50, 200 and 1000 tasks (4.03 s against 4.00 s at 1000 tasks and 2 MiB). At 8 MiB copying is 2 to 3 times slower (12.4 s against 4.7 s at 1000 tasks).

## Why per-element tasks always share

`mirai_map` queues every task up front, so a copied argument is held once per queued task. `loo_moment_match` (400 tasks, about 200 KB of arguments) used 10.5% more peak memory with copies (+117 MiB, interval +113.7 to +122.3 MiB) and the same time. Sharing every argument in that mode brought memory back to +0.2% (interval -0.1 to +0.4%). Tasks sent that way are expensive, so the shared-read cost does not matter.

## Checks

foil, 4 workers, 10 repetitions, always-share commit `b131a115` against the rule above: `diamonds_function` -21% time, `loo_function` level, `moment_matching` and `vig_moment_matching` level on time and memory, `psis_matrix_large`, `psis_matrix_small`, `relative_eff` and `subsample_methods` level.

The cutoff is from these two sweeps on one machine. It is not tuned to the benchmark suite.
