# dev-notes

Notes and benchmark results for the parallel (mirai + mori) work on loo. Not part of the package;
delete before the PR is merged.

| Path | What it is |
|---|---|
| `machine_specs.md` | Hardware and package versions the benchmarks ran on |
| `mizu-evaluation.md` | Why a mizu backend was tried and dropped, with all numbers |
| `mori-sharing.md` | The measurements behind when `loo_map` shares arguments with mori |
| `benchmark-checks.md` | WAIC old vs new, batching variants, repetition counts, sleep, data size |
| `foil_results/parallel-vs-baseline/` | foil: pre-parallel loo (`af50d00b`) vs the parallel candidate |
| `foil_results/mizu/` | mirai vs mizu: foil runs and `loo_map` microbenchmarks |

## foil_results
Each benchmark folder holds `report.txt` (baseline, candidate, change and credible intervals)
and `measurements.csv` (per-run times and memory).

- `parallel-vs-baseline/windows_w4/`: 9 benchmarks, Windows, 4 workers. `configs/` has the foil configs. New runs come from `bench/run.sh` (see `bench/README.md`) and land in `foil_results/<platform>/w<N>/`.
- `windows/w2`, `w4`, `w8`: all 26 benchmarks at 2, 4 and 8 workers, baseline `loo-v3.0.0` (`af50d00b`) against candidate `24c52ca5`, Windows, 20 repetitions each (10 for the seven slowest), with `SUMMARY.txt` and `environment.txt`.
  `windows/mori/w4`: six benchmarks at 4 workers comparing `loo_map()` sharing every argument with mori against copying them. The runs come from `bench/run.sh` (`bench/foil.toml`).
- `linux/w2`, `w4`, `w8`, `w24`, `w32`: the same 26 benchmarks on a TIDE node (one worker count per node), time only. `linux/mori/w4`: the six mori benchmarks. `linux/serial1/w1`: one mirai daemon against `lapply`. `linux/cluster-*.txt` describe the nodes.
- `plot_scaling.R relative|absolute` (in `foil_results`): scaling plots per platform, `scaling_{linux,windows}_{relative,absolute}.png`, benchmarks ordered from longest to shortest sequential time.
- `mizu/windows_w4/`: 19 benchmarks, mirai baseline `45479fb6` vs mizu candidate, Windows, 4 workers,
  with `SUMMARY.txt`. `quick/` has two quick runs. `microbench/` has `loo_map` timings on Windows and WSL.

`foil_results/RESULTS.md` collects every table.
