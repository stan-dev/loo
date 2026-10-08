# loo v3.0.0: documentation overview

> Last update: 2026-09-29.
>
> Author: Florence Bockting

Related notes: [loo-v3-cleanup-todos.md](loo-v3-cleanup-todos.md),
[loo-v3-related-prs-issues.md](loo-v3-related-prs-issues.md).

## 1. Tutorials (existing) — users

| Document | Name | Status | Branch / PR |
| :-- | :-- | :-- | :-- |
| `vignettes/loo2-example.Rmd` | Using the loo package | changed: `model_compare()` | `integrate-loo_compare` (#380) |
| `vignettes/loo2-with-rstan.Rmd` | Writing Stan programs for use with the loo package | changed: `model_compare()` | `integrate-loo_compare` (#380) |
| `vignettes/loo2-elpd.Rmd` | Holdout validation and K-fold cross-validation | changed: prose; `elpd()` calls left | `integrate-loo_compare` (#380) |
| `vignettes/loo2-large-data.Rmd` | Using leave-one-out cross-validation for large data | changed: 3 of 4 `loo_compare()` calls replaced | `integrate-loo_compare` (#380) |
| `vignettes/loo2-weights.Rmd` | Bayesian stacking and pseudo-BMA weights | changed: parallelism | `parallelization` (#378) |
| `vignettes/loo2-non-factorized.Rmd` | Leave-one-out cross-validation for non-factorized models | unchanged | `master` |
| `vignettes/loo2-lfo.Rmd` | Approximate leave-future-out cross-validation | unchanged | `master` |
| `vignettes/loo2-moment-matching.Rmd` | Avoiding model refits with moment matching | unchanged | `master` |
| `vignettes/loo2-mixis.Rmd` | Mixture IS leave-one-out cross-validation | unchanged | `master` |
| external: `https://users.aalto.fi/~ave/CV-FAQ.html` | Cross-validation FAQ | unchanged | `master` (`_pkgdown.yml`) |

## 2. Tutorials (new in loo v3) — users

| Document | Name | Status | Branch / PR |
| :-- | :-- | :-- | :-- |
| `vignettes/articles-online-only/overview-measures.Rmd` | Predictive model performance: Overview of measures | exists; placeholders left | `pred_measure` (#363) |
| `vignettes/articles-online-only/pred-measure-workflow.Rmd` | Predictive model performance: Predictive schemes with `pred_measure()` | exists; placeholder left | `pred_measure` (#363) |
| `vignettes/articles-online-only/model-comparison.Rmd` | Model comparison: Explanation of `model_compare()` | exists; placeholders left | `integrate-loo_compare` (#380) |
| not yet created | Model comparison: Case Study with `model_compare()` | planned | issue #401 |

## 3. Other — users

| Document | Name | Status | Branch / PR |
| :-- | :-- | :-- | :-- |
| `R/loo-glossary.R` | The `loo` glossary | extended; general terms planned | `integrate-loo_compare` (#380) |
| ArXiV preprint | Predictive measures: A formula reference | in progress; no arXiv ID yet | in Overleaf |

## 4. Project management — users and developers

| Document | Name | Status | Branch / PR |
| :-- | :-- | :-- | :-- |
| `vignettes/migration-guide.Rmd` | Migration guide | exists; four versions on four branches | `loo-v3.0.0` (#379), extended on #363, #380, #378 |
| `vignettes/articles-online-only/roadmap.Rmd` | Development roadmap | draft | `roadmap` (#396) |

