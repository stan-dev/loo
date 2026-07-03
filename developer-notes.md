# Developer Notes: `pred_measure` Feature

> **Status:** In Progress  
> **Base branch:** `loo-v3.0.0`  
> **Compare branch:** `add-pred-measure`  
> **Related PR:** [#363](https://github.com/stan-dev/loo/pull/363)  
> **Contributors:** @florence-bockting, @avehtari, @VisruthSK, @jgabry  
> **Last updated:** 2026-07-03

These notes document internal design decisions and ongoing work for the
`pred_measure` feature. This PR **adds** the new API on top of `loo-v3.0.0`;
it does not include other planned 3.0.0 changes (e.g. broader `compare()` /
`psislw()` removals) documented on the `loo-v3.0.0` branch.

For the merge summary, see the PR description ([`notes/pr-pred_measure.md`](notes/pr-pred_measure.md)).

---

## Status at a Glance

| Area | Status |
|------|--------|
| `pred_measure` API | Done (initial release) |
| `measure_*()` built-ins | Done |
| Scoring rules (`measure_rps`) | Done |
| Documentation | In progress |
| `group_ids` grouping | Not started |
| `loo_compare` integration | Not started |

---

## Scope of this PR (`add-pred-measure` vs `loo-v3.0.0`)

### Added (did not exist on `loo-v3.0.0`)

- `R/pred_measure.R` — `insample_pred_measure()`, `loo_pred_measure()`,
  `kfold_pred_measure()`, `test_pred_measure()`, `pred_measure()`
- `R/pred_measure-compute.R`, `R/pred_measure-helpers.R`,
  `R/pred_measure-builtin.R` — orchestration and `measure_*()` implementations
- `supported_measures_list()`, `ptw_log_pred_density()`
- S3 print methods for `pred_measure`, `loo_pred_measure`, `kfold_pred_measure`
- `vignettes/migration-guide.Rmd`
- Website-only articles: `overview-measures.Rmd`, `pred-measure-workflow.Rmd`
- Test suite + pre-fitted fixtures + `test_data_generation.R`

### Changed on existing code (implementations retained)

- `elpd()`, `crps()`, `scrps()`, `loo_crps()`, `loo_scrps()`,
  `loo_predictive_metric()` — deprecated with migration docs; **same APIs and
  implementations** (e.g. `crps(x, x2, y)` permutation estimator unchanged)
- `elpd()` — refactored to `.elpd_matrix_impl()` to avoid double deprecation
  warnings
- Minor doc cross-references in `compare.R`, `psislw.R`
- `NEWS.md`, `NAMESPACE`, `_pkgdown.yml`, pkgdown CI workflow

### Not in scope of this PR

- Removing deprecated `crps()` / `elpd()` functions
- Broader `loo-v3.0.0` refactoring (see `dev/migration-guide-loo-v3.Rmd` on
  that branch)

---

## Design decisions (resolved)

### D2: Scaled RPS — separate function vs. argument

**Decision:** Both. `measure_rps(..., scaled = TRUE)` is the implementation;
`measure_srps()` is a convenience wrapper. An explicit `srps` built-in name is
required so users can request both in one call:

```r
pred_measure(y = y, ypred = ypred, measure = c("rps", "srps"), ...)
```

Without a separate `srps` name, duplicate measure names would require custom
functions.

### Unified scoring rule in the new API

`measure_rps()` subsumes CRPS, RPS, SCRPS, and SRPS for the **new** workflow:

- Single `ypred` matrix (not two independent draw matrices)
- PWM/ECDF estimator (Aki's implementation; verified against Seth's derivations)
- Works for continuous and ordered categorical outcomes
- `measure_srps()` = scaled variant

Deprecated `crps()` / `scrps()` remain for backward compatibility with the
`x`, `x2` permutation-based API. They are not removed in this PR.

### ELPD and `ic` in the new API

Design choices **internal to `pred_measure`** (not a migration from
`loo-v3.0.0`, which has no `pred_measure`):

- When `ylp` is supplied, `elpd` is computed as the base summary
- `ic` is **not** included automatically; request via `measure = "ic"`
- `measure_elpd()` is separate from deprecated `elpd()` (different return type:
  `"measure"` vs `"elpd_generic"`)

### K-fold + categorical / multinomial models (brms)

- [x] Upstream fix for 3D `brms::kfold_predict()` output
  ([brms#1889](https://github.com/paul-buerkner/brms/issues/1889),
  fixed in [brms#1890](https://github.com/paul-buerkner/brms/pull/1890),
  merged to `main`)
- [x] Vignettes and pkgdown CI install brms from GitHub `master`
  (`vignettes/children/LOAD-BRMS-GITHUB.txt`, `.github/workflows/pkgdown.yaml`,
  `tests/testthat/data-for-tests/test_data_generation.R`)
- [ ] Drop GitHub brms pin once a CRAN release includes the #1890 fix
- [ ] Verify `kfold_pred_measure()` with categorical/multinomial examples
  end-to-end (penguins fixture exists; confirm test/doc coverage)

---

## Open decisions

### D1: Sign convention for pointwise estimates

- **Context:** Measures differ in orientation (`rps`: lower is better;
  `srps`: higher is better). Aligning orientations may help comparisons.
- **Options:** `lower_is_better`, `orientation = "utility" / "loss"`,
  `revert_sign` (currently internal on some `measure_*()` functions)
- **Decision:** *pending*

### D3: Handling of `r_eff`

- **Context:** [stan-dev/posterior#446](https://github.com/stan-dev/posterior/issues/446)
- **Question:** How should `r_eff` be handled in `pred_measure` workflows?
- **Decision:** *pending*

---

## Tasks

### Refactoring (within new API)

- [x] Add `measure_elpd()` for the new API
- [x] Refactor deprecated `elpd()` to `.elpd_matrix_impl()` (no double warnings)
- [x] Add `measure_rps()` as unified scoring-rule implementation (single `ypred`)
- [x] Keep deprecated `crps()` / `scrps()` with `x`, `x2` API (not removed)
- [x] In `*_pred_measure()`, compute `elpd` as base when `ylp` supplied;
  require explicit `measure = "ic"` for information criterion
- [ ] Provide an interface to `loo_compare` and verify consistency
- [ ] Resolve `r_eff` handling *(see D3)*

### Implementation

- [x] Consolidate CRPS/RPS/SCRPS/SRPS in `measure_rps()` for new API
- [x] Resolve `srps` naming *(see D2)*
- [x] brms k-fold categorical support unblocked upstream *(see above)*
- [ ] Drop GitHub brms pin after CRAN release

### Documentation

- [x] Migration guide (`vignettes/migration-guide.Rmd`)
- [x] Overview of measures (`overview-measures.Rmd`)
- [x] Workflow article (`pred-measure-workflow.Rmd`)
- [x] Online-only articles published via `_pkgdown.yml`
- [ ] Formula derivations article (`pred_measure-formulas.Rmd`)
- [ ] Detailed per-measure descriptions (derivations where appropriate)
- [ ] Extend glossary (`R/loo-glossary.R`) — measure, metric, score, utility, loss

### Grouping via `group_ids`

- [ ] Work out grouping scenarios for LOO and k-fold
- [ ] Implement in `kfold_pred_measure()`, `loo_pred_measure()`, `pred_measure()`

  > **Blocker:** Implementation approach uncertain; grouping scenarios must be
  > understood before coding begins.

---

## General questions

- Rename `ic` → `information_criteria` for clarity?
- Should `measure_elpd()` also return `ic`, or keep them separate?
- What defines class `"loo"` on measure objects? (e.g. deprecated `elpd_generic`
  inherits `"loo"`)
- Should `elpd` always be computed when `ylp` is supplied, or allow
  `loo_pred_measure()` for non-ELPD measures only?

---

## References & resources

| Resource | Link / contact |
|----------|----------------|
| Cross-validation FAQ | <https://users.aalto.fi/~ave/CV-FAQ.html> |
| `r_eff` in posterior (Aki) | <https://github.com/stan-dev/posterior/issues/446> |
| brms k-fold 3D fix | [brms#1889](https://github.com/paul-buerkner/brms/issues/1889), [brms#1890](https://github.com/paul-buerkner/brms/pull/1890) |
| `rps` derivation (Seth) | @florence-bockting, @avehtari |
| `rps` implementation (Aki) | @florence-bockting, @avehtari |
| K-fold vignette | <https://mc-stan.org/loo/articles/loo2-elpd.html> |
| Broader loo 3.0.0 migration | `loo-v3.0.0` branch, `dev/migration-guide-loo-v3.Rmd` |
