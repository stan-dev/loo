# Developer Notes: `loo` Refactoring & `pred_measure` Feature

> **Status:** In Progress  
> **Branch:** `pred_measure`  
> **Related PR / Issue:** `https://github.com/stan-dev/loo/pull/363`  
> **Contributors:** @florence-bockting, @avehtari, @VisruthSK, @jgabry
> **Last updated:** 2026-06-03

These notes document the internal discussions and ongoing progress on the `loo`
refactoring and the addition of `pred_measure` features. They serve as a
transparent, shared reference for all developers involved, covering the current
status, open tasks, resolved items, and outstanding decisions.

## Status at a Glance

| Area | Status |
|---|---|
| Refactoring | In Progress |
| Scores & Metrics Implementation | In Progress |
| Documentation | In Progress |
| Grouping via `group_ids` | Not Started |

## Open Decisions

*These questions need to be resolved before or alongside implementation.
Discuss, record the agreed outcome, and update the status below.*

### D1: Sign convention for pointwise estimates

+ **Context:** Some measures have opposite orientations (e.g., `rps`: lower is
better; `srps`: higher is better). It may be useful to align measures of
different orientations so that comparisons are meaningful.
+ **Question:** Should we allow the sign of pointwise estimates to be reversed
via a user-facing argument?
+ **Options under consideration:**
  - `lower_is_better = TRUE / FALSE`
  - `orientation = "utility" / "loss"`
  - `revert_sign = TRUE / FALSE`
+ **Decision:** *pending*

### D2: Scaled `rps`: separate function vs. argument
+ **Context:** We consider to use `rps` as a unified score for continuous and
ordered categorical inputs, subsuming `crps`, `scrps`, `srps`, and `rps`. 
As they can use the same underlying implementation (see implementation shared 
by Aki; contact Aki or Flo)
+ **Question:** Do we need an explicit `srps` function, or is a `scale = TRUE /
FALSE` argument inside `rps` sufficient?
+ **Decision:** Yes. We need `srps` for easy user-input of measures in `pred_measure`. Example: 
```
pred_measure(
  y = y, ypred = ypred, measure = c("rps", "srps"), ...
)
```
Without having the scaled version as explicit form, we would have to create a custom measure with a different name otherwise we would have duplicate names.

### D3: Handling of `r_eff`
+ **Context:** There is an open upstream discussion on `r_eff` handling in
`posterior` by Aki:
[stan-dev/posterior#446](https://github.com/stan-dev/posterior/issues/446).
+ **Question:** How should `r_eff` be handled in this package?
+ **Decision:** *pending*



## Tasks

### Refactoring

- [x] Remove the "old" `crps` implementation that requires two independent
  `ypred` inputs; retain only the new ECDF-based implementation.
- [ ] Align the currently two existing `elpd` implementations 
(old and new implementation).
- [x] Move `ic` out of "base measures"; it should not be included
  automatically but must be explicitly requested by the user.
- [ ] Provide an interface to `loo_compare` and verify consistency with it.
- [ ] Resolve handling of `r_eff` *(see Decision D3 above)*.

### Implementation of Scores & Metrics

- [x] Consolidate `crps`, `rps`, `scrps`, and `srps` into a unified `rps`
  function usable for continuous and ordered categorical inputs (including
  scaled variants).
  - [x] Use the implementation shared by Aki with Flo.
  - [x] Verify consistency with Seth's derivations.
  - [x] Resolve the `srps` question *(see Decision D2 above)*.
- [ ] Address the `kfold_pred_measure()` limitation with families that return
  3D objects (e.g., categorical models produce an `S × N × K` array from
  `brms::kfold_predict()`):
  - [x] Open a corresponding issue in `brms`: [Issue#1889](https://github.com/paul-buerkner/brms/issues/1889)
    + Created also a corresponding [PR#1890](https://github.com/paul-buerkner/brms/pull/1890)
  - [ ] Until the `brms` issue is resolved: document which families are not
    yet supported and return a `"not implemented"` error with a pointer to
    the `brms` issue.

### Documentation

- [ ] Write detailed descriptions (including derivations where appropriate)
  for individual measures.
- [ ] Write overview document with short descriptions of all measures >
  [`vignettes/articles-online-only/pred_measure-formulas.Rmd`]
- [ ] Write overview document explaining the new `pred_measure` family >
  [`vignettes/articles-online-only/pred_measure-family.Rmd`]
- [ ] Ensure overview documents are published as online-only articles.
- [ ] Extend glossary with key terminology (measure, metric, score, utility,
  loss, ...) > [`R/loo-glossary.R`]

### Grouping via `group_ids`

- [ ] Work out the different grouping scenarios for `loo` and `kfold`.
- [ ] Add `group_ids` support to `kfold_pred_measure()`,
  `loo_pred_measure()`, and `pred_measure()`.

  > **Blocker:** The implementation approach is still uncertain. The
  > grouping scenarios must be fully understood before any implementation
  > begins.

## General questions

+ the function name `ic` for the information criteria is very short and ambiguous;
shall we change to `information_criteria` as a more descriptive name?
+ do we want to maintain a "generic elpd" measure that returns elpd and ic 
information? 
+ What is the criteria whether something is of class loo? For example, elpd.generic
is of class loo, why?

## References & Resources

| Resource | Link / Contact |
|---|---|
| Vignette: `kfold` and test data with `loo` | <https://mc-stan.org/loo/articles/loo2-elpd.html> |
| Discussion: `r_eff` handling in `posterior` (Aki) | <https://github.com/stan-dev/posterior/issues/446> |
| `rps` derivation by Seth | Contact @florence-bockting or @avehtari |
| `rps` implementation by Aki | Contact @florence-bockting or @avehtari |
