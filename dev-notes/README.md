# dev-notes

Notes and benchmark results for the parallel (mirai + mori) work on loo. Not part of the package;
delete before the PR is merged.

| Path | What it is |
|---|---|
| `machine_specs.md` | Hardware and package versions the benchmarks ran on |
| `mizu-evaluation.md` | Why a mizu backend was tried and dropped, with all numbers |
| `foil_results/parallel-vs-baseline/` | foil: pre-parallel loo (`af50d00b`) vs the parallel candidate |
| `foil_results/mizu/` | mirai vs mizu: foil runs and `loo_map` microbenchmarks |

## foil_results
Each benchmark folder holds `report.txt` (baseline, candidate, change and credible intervals)
and `measurements.csv` (per-run times and memory).

- `parallel-vs-baseline/windows_w4/`: 9 benchmarks, Windows, 4 workers. `configs/` has the foil configs.
- `mizu/windows_w4/`: 19 benchmarks, mirai baseline `45479fb6` vs mizu candidate, Windows, 4 workers,
  with `SUMMARY.txt`. `quick/` has two quick runs. `microbench/` has `loo_map` timings on Windows and WSL.

Only these runs exist; earlier 2- and 8-worker and WSL foil runs were not kept.
