# mizu as a loo_map backend: evaluation

Question: does running `loo_map()` on a local mizu pool (shared-memory, work-stealing) beat the
existing mirai + mori path enough to justify supporting it?

Verdict: no. The mizu backend was implemented, tested, benchmarked and reverted. The code is in git history
(commits `cbcff0cb`, `c4bb23eb`, `93386aaa`; reverted in `e824b2ed`) on branch `parallel`.

## What was built
- `loo_map()` used `mizu::mizu_map()` when `mizu::mizu_default_pool()` was set, else mirai (unchanged).
- `batch = TRUE`: one block of consecutive elements per worker. `batch = FALSE`: one element per task.
  Arguments went through `mori::share()`. Worker warnings were captured per block and re-signalled on the host.
- Full test suite passed on both backends (`LOO_TEST_BACKEND=mizu` switched the suite onto a mizu pool).
  Results from mizu, mirai and sequential runs were identical.

## Why it does not help loo
- mizu's published speedups (13x data frame round trip, 37x over mirai on 5,000 ~10 us tasks) are against
  mirai sending one task per element and serializing large objects over sockets.
  loo's baseline already sends one block per daemon and shares data with mori, so that overhead is gone.
- In `loo(llfun, ...)` with 4000 observations, 98-99% of wall time is inside `loo_map()` on both backends, at
  about 2 ms of real work per element. Dispatch is a few milliseconds of a multi-second run.
- Per-element warning capture is needed because worker warnings do not propagate. Doing it per element costs
  more than it saves for cheap elements (100k cheap elements: 315 ms per element vs 199 ms in blocks),
  `globalCallingHandlers()` cannot be installed inside a task. Hence the blocks.

## Results
Raw files are under `foil_results/mizu/`.

### loo_map microbenchmark, median ms, 4 workers (`foil_results/mizu/microbench/`)
| case | Windows mirai | Windows mizu | WSL mirai | WSL mizu |
|---|---|---|---|---|
| call_8_cheap (small repeated call) | 0.53 | 0.43 | 3.09 | 0.59 |
| tiny_2000 | 3.23 | 3.35 | 7.78 | 4.27 |
| psislike_4000 (compute-bound) | 306 | 347 | 469 | 343 |
| cheap_100k | 220 | 240 | 232 | 355 |
| skewed_1000 (per element) | 237 | 160 | 585 | 119 |
| per_element_200 (per element) | 37 | 14 | 139 | 24 |

Windows 2 and 8 workers show the same pattern (per element 1.3-2.7x faster on mizu, compute-bound tied).
Memory (private bytes summed over host and workers) was comparable.
On WSL, mirai crashed with "Connection reset" in the 1000-task per-element case when run together with
the other cases; mizu did not.

### End to end (`microbench/wsl_loo_function_end_to_end.txt`)
Tie within noise on Windows and WSL.

### Foil, mirai baseline `45479fb6` vs mizu candidate, Windows, 4 workers (`foil_results/mizu/windows_w4/`, 19 benchmarks)
Most benchmarks within +-4% (`foil_results/mizu/windows_w4/SUMMARY.txt`). Faster on mizu: wells_loo_subsample -5.6%,
vig_moment_matching -3.8%, diamonds_matrix -1.9%. Slower: radon_loo_function +6.6%, vig_milk_kline +4.0%,
brms_toplevel +3.8%, loo_five_calls +3.7%. mizu's pool start (about 0.2 s, blocking) counts in every whole-process time.
Quick runs (`foil_results/mizu/quick/`): moment_matching -16% in the first quick run, loo_function +4.5%.
2 and 8 worker foil runs were not completed.

## Where mizu would pay off
Only when each element costs under about 100 us, or for many tiny repeated calls on Linux. loo does not
produce that regime. Making loo faster means cheaper per-observation PSIS / ESS work, not a different backend.

