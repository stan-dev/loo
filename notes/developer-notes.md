Developer Notes: `pred_measure` Feature
================

- [Status at a Glance](#status-at-a-glance)
- [Scope of this PR (`add-pred-measure` vs
  `loo-v3.0.0`)](#scope-of-this-pr-add-pred-measure-vs-loo-v300)
- [Design decisions (resolved)](#design-decisions-resolved)
- [Open decisions](#open-decisions)
- [Tasks](#tasks)
- [General questions](#general-questions)
- [References & resources](#references--resources)
- [Appendix: Numerical comparisons (deprecated vs new
  API)](#appendix-numerical-comparisons-deprecated-vs-new-api)

> **Status:** In Progress  
> **Base branch:** `loo-v3.0.0`  
> **Compare branch:** `pred_measure` (+ `integrate-loo_compare`)  
> **Related PR:** [\#363](https://github.com/stan-dev/loo/pull/363)  
> **Contributors:** @florence-bockting, @avehtari, @VisruthSK, @jgabry  
> **Last updated:** 2026-07-07

These notes document internal design decisions and ongoing work for the
`pred_measure` feature. This PR **adds** the new API on top of
`loo-v3.0.0`; it does not include other planned 3.0.0 changes
(e.g. broader `compare()` / `psislw()` removals) documented on the
`loo-v3.0.0` branch.

For the merge summary, see the PR description
([`internal-notes/pr-pred_measure.md`](../internal-notes/pr-pred_measure.md)).

------------------------------------------------------------------------

## Status at a Glance

| Area                          | Status                 |
|-------------------------------|------------------------|
| `pred_measure` API            | Done (initial release) |
| `measure_*()` built-ins       | Done                   |
| Scoring rules (`measure_rps`) | Done                   |
| Documentation                 | In progress            |
| `group_ids` grouping          | Not started            |
| `loo_compare` integration     | Done (`integrate-loo_compare`) |

------------------------------------------------------------------------

## Scope of this PR (`add-pred-measure` vs `loo-v3.0.0`)

### Added (did not exist on `loo-v3.0.0`)

- `R/pred_measure.R` — `insample_pred_measure()`, `loo_pred_measure()`,
  `kfold_pred_measure()`, `test_pred_measure()`, `pred_measure()`
- `R/pred_measure-compute.R`, `R/pred_measure-helpers.R`,
  `R/pred_measure-builtin.R` — orchestration and `measure_*()`
  implementations
- `supported_measures_list()`, `ptw_log_pred_density()`
- S3 print methods for `pred_measure`, `loo_pred_measure`,
  `kfold_pred_measure`
- `vignettes/migration-guide.Rmd`
- Website-only articles: `overview-measures.Rmd`,
  `pred-measure-workflow.Rmd`
- Test suite + pre-fitted fixtures + `test_data_generation.R`
- `loo_compare()` multi-measure path for `loo_pred_measure` objects
  (`integrate-loo_compare`)

### Changed on existing code (implementations retained)

- `elpd()`, `crps()`, `scrps()`, `loo_crps()`, `loo_scrps()`,
  `loo_predictive_metric()` — deprecated with migration docs; **same
  APIs and implementations** (e.g. `crps(x, x2, y)` permutation
  estimator unchanged)
- `elpd()` — refactored to `.elpd_matrix_impl()` to avoid double
  deprecation warnings
- Minor doc cross-references in `compare.R`, `psislw.R`
- `loo_compare()` — extended for `loo_pred_measure` objects: `rank_by`,
  multi-measure paired diffs, updated `print.compare.loo(measures = ...)`;
  classic `loo` path unchanged
- `R/loo-glossary.R` — multi-measure comparison columns (`{measure}_diff`,
  `rank_by`, etc.)
- `NEWS.md`, `NAMESPACE`, `_pkgdown.yml`, pkgdown CI workflow

------------------------------------------------------------------------

## Design decisions (resolved)

### D2: Scaled RPS — separate function vs. argument

**Decision:** Both. `measure_rps(..., scaled = TRUE)` is the
implementation; `measure_srps()` is a convenience wrapper. An explicit
`srps` built-in name is required so users can request both in one call:

``` r
pred_measure(y = y, ypred = ypred, measure = c("rps", "srps"), ...)
```

Without a separate `srps` name, duplicate measure names would require
custom functions.

### Unified scoring rule in the new API

`measure_rps()` subsumes CRPS, RPS, SCRPS, and SRPS for the **new**
workflow:

- Single `ypred` matrix (not two independent draw matrices)
- PWM/ECDF estimator (Aki’s implementation; verified against Seth’s
  derivations)
- Works for continuous and ordered categorical outcomes
- `measure_srps()` = scaled variant

Deprecated `crps()` / `scrps()` remain for backward compatibility with
the `x`, `x2` permutation-based API.

### ELPD and `ic` in the new API

Design choices **internal to `pred_measure`** (not a migration from
`loo-v3.0.0`, which has no `pred_measure`):

- When `ylp` is supplied, `elpd` is computed as the base summary
- `ic` is **not** included automatically; request via `measure = "ic"`
- `measure_elpd()` is separate from deprecated `elpd()` (different
  return type: `"measure"` vs `"elpd_generic"`)

### K-fold + categorical / multinomial models (brms)

- [x] Upstream fix for 3D `brms::kfold_predict()` output
  ([brms#1889](https://github.com/paul-buerkner/brms/issues/1889), fixed
  in [brms#1890](https://github.com/paul-buerkner/brms/pull/1890),
  merged to `main`)
- [x] Vignettes and pkgdown CI install brms from GitHub `master`
  (`vignettes/children/LOAD-BRMS-GITHUB.txt`,
  `.github/workflows/pkgdown.yaml`,
  `tests/testthat/data-for-tests/test_data_generation.R`)
- [ ] Drop GitHub brms pin once a CRAN release includes the \#1890 fix
- [ ] Verify `kfold_pred_measure()` with categorical/multinomial
  examples end-to-end (penguins fixture exists; confirm test/doc
  coverage)

### D4: `loo_compare()` for `loo_pred_measure` objects

**Decision:** Extend existing `loo_compare()`

- When all inputs are `loo_pred_measure` objects, compute paired
  differences for every measure common to all models
- Rank models by `rank_by` (default `"elpd"`); top-ranked model is the
  reference for all `{measure}_diff` columns
- ELPD-family measures keep `elpd_diff` / `se_diff`; other measures use
  `{measure}_diff` / `{measure}_se_diff`
- `p_worse` and `diag_diff` apply to ELPD only; `diag_elpd` per model as
  before
- Loss measures (MSE, RMSE, MAE, IC, Brier score, RPS) compared on a common
  utility scale (higher is better): sign flipped from the raw loss orientation
  so worse models have negative diffs, consistent with ELPD. Orientation is
  read from the `loss` element of `measure_info` on each
  `*_pred_measure()` result; attribute `sign_converted_measures` records
  affected measures. A short message is emitted at compare time; full
  interpretation is in `?loo_compare` / `?loo-glossary`.
- Pointwise SEs use the same paired formula as ELPD when the overall
  estimate is a sum or mean of pointwise contributions; otherwise
  `{measure}_se_diff` is `NA` (e.g. `r2`, `mse`, `rmse`)
- Reuse `elpd_diffs`, `se_elpd_diff`, `diag_diff`, `diag_elpd`, and
  many-model order-statistic check (with `rank_by` when applicable)
- `print.compare.loo(measures = ...)` shows one or all measure diff tables

Implemented on branch `integrate-loo_compare`; tests in `test_compare.R`
with fixture `test_data_roaches_compare.Rds`.

------------------------------------------------------------------------

## Open decisions

### D1: Sign convention for pointwise estimates

- **Context:** Measures differ in orientation (e.g. ELPD/CRPS on a utility
  scale; MSE and Brier score as losses). `loo_compare()` aligns them for
  paired differences.
- **Decision:** Values are always stored on the measure's own scale. A single
  `loss` flag records the orientation --- `.measure_spec` for built-ins,
  `attr(fun, "measure_loss")` for custom measures --- and is recorded per
  measure in the `measure_info` attribute of each `*_pred_measure()` result.
  `model_compare()` sign-flips measures with `loss = TRUE` for utility-scale
  `{measure}_diff`. The per-call `higher_is_better` control has been removed:
  it only changed the sign of stored values, and comparisons were invariant
  to it.

### D3: Handling of `r_eff`

- **Context:**
  [stan-dev/posterior#446](https://github.com/stan-dev/posterior/issues/446)
- **Question:** How should `r_eff` be handled in `pred_measure`
  workflows?
- **Decision:** *pending*

### D5: Weighted `E|X − X'|` — why the derivation in `crps_pwm.pdf` is not used

- **Decision:** *resolved.* `.exx_pwm()` (`R/pred_measure-helpers.R`) estimates
  `E|X − X'|` with the bias-corrected weighted Gini mean difference

  ``` r
  EXX = sum_{i != j} w_i w_j |x_i - x_j| / (1 - sum_i w_i^2)
  ```

  computed on sorted draws as `2 * sum(w_s * x_s * (2 * C_s - w_s - 1)) / (1 - sum(w^2))`,
  with `C_s = sum_{k <= s} w_k`. At equal weights this is exactly the unbiased
  PWM estimator of Taillardat et al. (2016) with the `1 / (S (S - 1))`
  normalization, so the weighted and unweighted paths agree.

- **Why not `notes/crps_pwm.pdf` (section 0.3).** The note derives a weighted
  estimator from `E|X − X'| = 2 (E[X] − E[X_{1,1:2}])`, taking the probability
  that the `s`-th order statistic is the smaller of a random pair to be
  `2 * w_s * (1 - C_s) / (1 - w_s)`, where the factor 2 is said to handle order
  invariance. Under weighted sampling *without replacement* the two orderings do
  not have the same probability:

  ```
  P(x_s is the pair minimum) = w_s (1 - C_s) / (1 - w_s)        # x_s drawn first
                             + sum_{k > s} w_k w_s / (1 - w_k)  # x_s drawn second
  ```

  The two lines coincide only when all weights are equal — which is why the
  note's section 0.2 (unweighted) is exact and section 0.3 is not. Doubling the
  first line overweights draws that are both heavy and small, and the implied
  coefficients no longer sum to zero. `EXX` is then neither shift invariant nor
  guaranteed non-negative, and `measure_srps()` takes `log()` of a negative
  number. (Deriving the exact pair probabilities is not a small fix: they are
  the second-order inclusion probabilities of PPS sampling without replacement,
  which do not factorize in general.)

  Measured on a toy sample of `S = 10` standard normal draws with weights
  `(0.02, ..., 0.02, 0.82)`, and on the roaches fixture
  (`tests/testthat/data-for-tests/test_data_roaches.Rds`, 262 obs × 400 draws):

  | case | note (0.3) | implemented |
  |---|---|---|
  | toy `EXX` | 2.647 | 1.547 |
  | toy `EXX`, draws shifted by +100 | 136.525 | 1.547 |
  | roaches obs 230 (max weight 0.833) | −44.67 → `NaN` in `srps` | 6.53 |
  | roaches obs 16 (max weight 0.778) | 392.25 | 22.97 |

  The note's coefficients sum to 1.339 rather than 0 in the toy case. The
  unweighted path is unaffected — there the note is exact and agrees with the
  implementation to machine precision.

- **Independent backing (ArviZ).** `arviz-stats` implements the same weighted
  PWM score in the same PSIS-LOO setting, in
  [`_loo_score()`](https://github.com/arviz-devs/arviz-stats/blob/main/src/arviz_stats/base/diagnostics.py):

  ``` python
  f_minus = cumulative_weights - weights_sorted
  bracket = 2.0 * f_minus + weights_sorted - 1.0
  gini_mean_difference = 2.0 * np.sum(weights_sorted * values_sorted * bracket)
  ```

  This is our numerator exactly (`2 * C_s - w_s - 1 = 2 * f_minus + w_s - 1`),
  and it carries no `1 / (1 - w_s)` factor — i.e. ArviZ independently arrived at
  the weighted Gini mean difference rather than at the note's estimator. The
  same double-sum form is the standard survey-weighted Gini estimator,
  `sum_k sum_l w_k w_l |y_k - y_l| / (2 N̂ Ŷ)`.

- **Why we keep the `1 / (1 - sum w^2)` correction that ArviZ omits.** ArviZ
  computes the plug-in version; ours divides by `1 - sum_i w_i^2 = 1 - 1/S_eff`,
  the standard reliability-weights bias correction (the weighted-variance
  analogue), which at equal weights is the `(S - 1) / S` "fair score" correction
  of Ferro (2014) discussed by Zamo & Naveau (2018). Two reasons:

  1. **The bias is per observation, not a constant.** With equal weights the
     plug-in is low by a fixed `(S - 1) / S`, which cancels everywhere. With PSIS
     weights the factor is `1 - 1/S_eff`, and on the roaches fixture `S_eff`
     ranges from 400 (median 368, factor 0.9973) down to 1.4 (factor 0.2973).
     27 of 262 observations differ by more than 1%, and in `srps` the omitted
     correction lands as an additive per-observation shift in `-0.5 * log(EXX)`
     of up to 0.61 — largest exactly where the importance weights are already
     concentrated.
  2. **Consistency with the unweighted path.** `measure_rps()` uses the unbiased
     PWM estimator when `log_weights` is `NULL`. Without the correction, uniform
     `log_weights` would no longer reproduce that result (off by `(S - 1) / S`);
     the test *"uniform log-weights reproduce the unweighted measure_rps()"*
     asserts that they do.

  Adopting ArviZ's plug-in form would therefore mean changing the unweighted
  path as well, which changes published `measure_rps()` output and drops the
  fair-score correction that `crps()`'s own references argue for.

- **Follow-up:** report the section 0.3 issue to the author of
  `notes/crps_pwm.pdf`; the note's unweighted result stands, only the weighted
  generalization needs revising. Note also that the printed code for
  `EXX_compute_pwm()` in section 0.4 has a typo (`- 2` should be
  `- 2 * (S + 1) / (S - 1)`); it disagrees with the note's own formula and is
  not shift invariant. `.exx_pwm()` follows the formula, not that snippet.

------------------------------------------------------------------------

## Tasks

### Refactoring (within new API)

- [x] Add `measure_elpd()` for the new API
- [x] Refactor deprecated `elpd()` to `.elpd_matrix_impl()` (no double
  warnings)
- [x] Add `measure_rps()` as unified scoring-rule implementation (single
  `ypred`)
- [x] Keep deprecated `crps()` / `scrps()` with `x`, `x2` API unchanged
- [x] In `*_pred_measure()`, compute `elpd` as base when `ylp` supplied;
  require explicit `measure = "ic"` for information criterion
- [x] Document and test deprecated vs new API comparisons *(see
  appendix)*
- [x] Provide an interface to `loo_compare` and verify consistency
- [ ] Decide whether the `loo_compare` S3 tree stays. `loo_compare.default`
  (`R/loo_compare.R:34`) and `loo_compare.psis_loo_ss_list`
  (`R/loo_compare.R:46`) are now two-line pass-throughs to their
  `model_compare` counterparts, so the whole tree may be a thin back-compat
  shim. Either keep it deliberately, as with `old_nms` / `convert_old_object()`,
  or drop it as a set — not one method at a time.
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
- [x] Extend glossary (`R/loo-glossary.R`) — multi-measure `loo_compare` columns
- [ ] Extend glossary further — measure, metric, score, utility, loss
  (general terms)

### Grouping via `group_ids`

- [ ] Work out grouping scenarios for LOO and k-fold

- [ ] Implement in `kfold_pred_measure()`, `loo_pred_measure()`,
  `pred_measure()`

  > **Blocker:** Implementation approach uncertain; grouping scenarios
  > must be understood before coding begins.

------------------------------------------------------------------------

## General questions

- Rename `ic` → `information_criteria` for clarity?
- Should `measure_elpd()` also return `ic`, or keep them separate?
- What defines class `"loo"` on measure objects? `loo_pred_measure` inherits
  `"loo"` (see `integrate-loo_compare`); deprecated `elpd_generic` also
  inherits `"loo"`.
- Should `elpd` always be computed when `ylp` is supplied, or allow
  `loo_pred_measure()` for non-ELPD measures only?

------------------------------------------------------------------------

## References & resources

| Resource                    | Link / contact                                                                                                               |
|-----------------------------|------------------------------------------------------------------------------------------------------------------------------|
| Cross-validation FAQ        | <https://users.aalto.fi/~ave/CV-FAQ.html>                                                                                    |
| `r_eff` in posterior (Aki)  | <https://github.com/stan-dev/posterior/issues/446>                                                                           |
| brms k-fold 3D fix          | [brms#1889](https://github.com/paul-buerkner/brms/issues/1889), [brms#1890](https://github.com/paul-buerkner/brms/pull/1890) |
| `rps` derivation (Seth)     | @florence-bockting, @avehtari                                                                                                |
| `rps` implementation (Aki)  | @florence-bockting, @avehtari                                                                                                |
| K-fold vignette             | <https://mc-stan.org/loo/articles/loo2-elpd.html>                                                                            |
| Broader loo 3.0.0 migration | `loo-v3.0.0` branch, `dev/migration-guide-loo-v3.Rmd`                                                                        |

------------------------------------------------------------------------

## Appendix: Numerical comparisons (deprecated vs new API)

Simulations and tests that validate design choices in this PR. Automated
comparisons live in `tests/testthat/test_crps.R` (CRPS/RPS section).

### CRPS / RPS

Deprecated `crps()` / `scrps()` and new `measure_rps()` /
`measure_srps()` target the same scoring rules but use different
estimators.

#### API mapping

| Deprecated        | New workflow                                    | Notes                               |
|-------------------|-------------------------------------------------|-------------------------------------|
| `crps(x, x2, y)`  | `-measure_rps(y, ypred = x)$pointwise` | `crps()` returns the negated unscaled score |
| `scrps(x, x2, y)` | `measure_srps(y, ypred = x)`                    | Same sign convention                |
| `loo_crps(...)`   | `loo_pred_measure(..., measure = "rps")`        | Additional LOO weighting difference |
| `loo_scrps(...)`  | `loo_pred_measure(..., measure = "srps")`       | Additional LOO weighting difference |

#### Sources of numerical difference

1.  **Sign convention (unscaled only).** `measure_rps()` returns
    `EXy − 0.5·EXX`, the Gneiting & Raftery (2007) loss (lower is better);
    `crps()` returns its negation `0.5·EXX − EXy`. Negate `measure_rps()`
    to match `crps()`. Scaled scores (`scrps` / `measure_srps`) already
    share the formula `−EXy/EXX − 0.5·log(EXX)` and are utilities.

2.  **EXX estimator (in-sample).** Both estimate `E|X − X'|`, but:

    - **Deprecated:** two draw matrices `x`, `x2`; one random shuffle
      per permutation (`EXX_compute()` in `R/crps.R`).
    - **New:** single `ypred` matrix; PWM on sorted draws (Zamo &
      Naveau, 2018).

    `EXy = E|X − y|` is **identical** between APIs. After sign
    alignment, all pointwise differences come from EXX.

3.  **LOO weighting.** `loo_crps()` shuffles a second draw matrix and
    applies joint PSIS weights to `|x − x2|`.
    `loo_pred_measure(..., measure = "rps")` uses weighted PWM on a
    single `ypred` with PSIS weights from `ylp` only — so LOO
    differences combine EXX method and importance-weighting approach.

    For why the weighted PWM estimator does not follow section 0.3 of
    `notes/crps_pwm.pdf`, see D5 in *Open decisions*.

#### Key results (reference simulation)

<figure>
<img src="figures/crps-rps-comparison.png"
alt="CRPS/RPS comparison: PWM vs permutation EXX estimators" />
<figcaption aria-hidden="true">CRPS/RPS comparison: PWM vs permutation
EXX estimators</figcaption>
</figure>

*Figure: left — per-observation EXX estimates at the brms default draw
count (S = 4000 post-warmup draws: 4 chains × 1000); right — mean
relative EXX error decreases with more posterior draws.*

| Metric (S = 100, n = 30)                   |  Value |
|:-------------------------------------------|-------:|
| Mean rel. \|EXX_perm - EXX_pwm\| / EXX_pwm | 0.0688 |
| cor(pointwise CRPS, sign-aligned)          | 0.9923 |
| Mean \|pointwise CRPS diff\|               | 0.0339 |
| cor(pointwise SCRPS, SRPS)                 | 0.9911 |

#### In-sample outcome comparison

<figure>
<img src="figures/crps-rps-outcomes.png"
alt="CRPS/RPS outcome comparison across replications" />
<figcaption aria-hidden="true">CRPS/RPS outcome comparison across
replications</figcaption>
</figure>

*Figure: 200 simulations (S = 100, n = 30). Left — `crps()` vs
`-measure_rps()`; right — `scrps()` vs
`measure_srps()`.*

#### LOO outcome comparison

<figure>
<img src="figures/loo-crps-rps-outcomes.png"
alt="LOO CRPS/RPS outcome comparison across replications" />
<figcaption aria-hidden="true">LOO CRPS/RPS outcome comparison across
replications</figcaption>
</figure>

*Figure: 200 LOO simulations (S = 100, n = 30). Left — `loo_crps()` vs
sign-aligned `loo_pred_measure(..., measure = "rps")`; right —
`loo_scrps()` vs `loo_pred_measure(..., measure = "srps")`.*

**CRPS/RPS takeaway:** Results are highly correlated but not
interchangeable. The PWM estimator is lower-variance and requires only
one draw matrix. For migration, compare trends and rankings rather than
expecting pointwise equality.

### ELPD / IC

Deprecated `elpd()` and new `measure_elpd()` / `measure_ic()` compute
the same in-sample pointwise log predictive density and information
criterion from a log-likelihood matrix. The difference is return
structure (`elpd_generic` with both columns vs separate `"measure"`
objects), not the underlying formula.

| Deprecated                      | New workflow        | Notes                             |
|---------------------------------|---------------------|-----------------------------------|
| `elpd(ylp)$estimates["elpd", ]` | `measure_elpd(ylp)` | Same `lppd_i` computation         |
| `elpd(ylp)$estimates["ic", ]`   | `measure_ic(ylp)`   | `ic_i = -2 * lppd_i` in both APIs |

#### Outcome comparison

<figure>
<img src="figures/elpd-ic-outcomes.png"
alt="ELPD / IC outcome comparison across replications" />
<figcaption aria-hidden="true">ELPD / IC outcome comparison across
replications</figcaption>
</figure>

*Figure: 200 simulations (S = 100, n = 30). Left — `elpd()` vs
`measure_elpd()`; right — `ic` from `elpd()` vs `measure_ic()`.*

**ELPD/IC takeaway:** In-sample estimates match between deprecated and
new APIs. Migration is about return type and `*_pred_measure()` workflow
integration, not numerical differences.
