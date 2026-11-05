# loo v3.0.0: related PRs and issues

> Last update: 2026-09-29
>
> Author: Florence Bockting

This note lists the PRs and issues that relate to the v3.0.0 refactoring.
The open work is in [loo-v3-cleanup-todos.md](loo-v3-cleanup-todos.md).

## Merge chain

| PR | Head branch | Base branch | Brings |
| :-- | :-- | :-- | :-- |
| #380 | `integrate-loo_compare` | `pred_measure` | `model_compare()`, `model-comparison.Rmd`, `notes/design-discussions/`, compare fixture, deprecation of `loo_compare()` |
| #363 (draft) | `pred_measure` | `loo-v3.0.0` | `*_pred_measure()`, `measure_*()`, `overview-measures.Rmd`, `pred-measure-workflow.Rmd`, deprecations of `elpd()`, `crps()` and others, test fixtures |
| #379 (draft) | `loo-v3.0.0` | `master` | the release |
| this | `loo-v3-cleanup` | `loo-v3.0.0` | the items below |

**Targeted for v3.0.x:**

| PR | Head branch | Base branch | Brings |
| :-- | :-- | :-- | :-- |
| #378 (draft) | `parallelization` | `loo-v3.0.0` | mirai/mori parallelism, fixes #308. Not part of v3.0.0, but it merges into the same base. |

## Issues that the refactoring fixes

| Issue | Title | Fixed by |
| :-- | :-- | :-- |
| #281 | New functions for better support for different scores and metrics | #363 |
| #223 | (loo_)(s)crps could ask for only one argument with predictions | #363 |
| #213 | loo_predictive_metric and loo_crps could accept psis objects | #363 |
| #201 | Add R2 | #363 |
| #135 | user defined loss/utility functions | #363 |
| #220 | loo_compare for crps and loo_crps | #363 (part), #380 (rest) |

All six are open. GitHub closes an issue from a "Fixes #n" line only when the
PR merges into `master`. #363 merges into `loo-v3.0.0`.

## Open PRs into `master` that touch the same code

| PR | Topic | Relation to this refactoring |
| :-- | :-- | :-- |
| #398 | `crps.numeric()` passes `permutations` wrongly (fixes #397) | Edits `R/crps.R`. #363 deprecates `crps()` in the same file. Merge conflict. |
| #393 (draft) | `loo_compare()` returns a data.frame for subsampling (fixes #392) | #380 deprecates `loo_compare()`. |
| #399 (draft) | Correct `pred_measure` values after post-hoc methods | |
| #340 | Export `srs_diff_est()` (fixes #333) | Subsampling comparison. `model_compare()` has a subsampling method. |
| #178 (draft) | LOO difference plot | #394 (milestone v3.0.0) extends it to the new measures. |
| #291 (draft) | Format the project with Air | Rewrites every file. Merge it before #363 or after #379, never between them. |
| #396 | Roadmap on the website | Describes the v3.0.0 plan. Update it when this list changes. |

## Other issues, by target release

All of these issues have the GitHub milestone v3.0.0. This list gives the
release that each one targets.

**v3.0.0**

+ #401 (case study for model comparison)
+ #388 (clean up issues): general clean-up and release preparation
+ #353 (warnings and messages): standardize and improve warning/messaging behavior implemented in loo

**v3.0.0, developer discussion**

+ #384 (LLM code optimizations): How do we want to integrate this in the development process

**v3.0.0 and later**

+ #343 (subset psis objects)
+ #249 (more `posterior` functions): avoid duplicated code and decide which functions should go into posterior and which should stay in loo
+ #394 (extend the difference plot)

**v3.0.x**

+ #308 / #378 (parallelism): improve parallelization in loo

**v3.x.x**

+ (loo diagnostics): improved posthoc diagnostics using loo
+ #192 (n_eff to ESS): improve diagnostics

**Target not decided**

+ #385 (mutation testing)

