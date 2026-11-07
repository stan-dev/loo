# Vignette benchmarks

Times only the loo calls that vignettes and case studies run, not the model fits. The log-likelihood draws come from fits made once beforehand, so each workload is: `relative_eff()`, `loo(save_psis = TRUE)`, `waic()` where the vignette calls it, and `loo_compare()`.

| Workload | Source | Models | Draws x observations |
|---|---|---|---|
| `bayesian-workflow/roaches` | [Bayesian Workflow](https://avehtari.github.io/Bayesian-Workflow), `roaches/roaches.R` | Poisson, negative binomial, Poisson with varying intercepts | 4000 x 262 |
| `bayesian-workflow/arsenic_loo` | Bayesian Workflow, `misc/chapter_08/section_08_04/arsenic_loo.R` | two logistic regressions on the wells data | 4000 x 3020 |
| `loo2-example` | `vignettes/loo2-example.Rmd` | two roaches models | 4000 x 262 |
| `loo2-with-rstan` | `vignettes/loo2-with-rstan.Rmd` | two wells models | 4000 x 3020 |
| `loo2-weights` | `vignettes/loo2-weights.Rmd` | milk and Kline model sets | small |

The other loo vignettes run through foil in `bench/foil.toml`: `vig_large_data` is `loo2-large-data`, `vig_moment_matching` is `loo2-moment-matching`, `vig_lfo_psis_loop` is `loo2-lfo`. `loo2-elpd`, `loo2-mixis` and `loo2-non-factorized` are not benchmarked.

## Running
`fit.R <Bayesian-Workflow checkout> <output dir>` fits the book models with rstanarm and saves the log-likelihood arrays. `time.R <base|cand> <that output dir> <directory holding real_roach.rds, real_wells.rds, real_milk_kline.rds>` times one implementation: `base` is pre-parallel loo (`cores`), `cand` is this branch (`parallel`, daemons started once beforehand). Both read libraries from `~/lx/libb` and `~/lx/libc`.

`results.csv` holds the timings: median and 10th and 90th percentile of 7 repetitions, in seconds.
