# Design: post-hoc LOO corrections in `pred_measure`

| | |
|---|---|
| **Package / version** | `loo` 3.0.0; `brms` 2.23.2 |
| **Author(s)** | Florence Bockting |
| **Status** | Draft. |
| **Created / last updated** | 2026-09-25 |
| **Related** | branch `integrate-loo_compare`; vignette `pred-measure-workflow.Rmd` |

## Contents

1. [Summary](#1-summary)
2. [Background and motivation](#2-background-and-motivation)
3. [Goals and non-goals](#3-goals-and-non-goals)
4. [User-facing design (API)](#4-user-facing-design-api)
   - [4.1 Usage scenarios](#41-usage-scenarios)
   - [4.2 Proposed interface](#42-proposed-interface)
   - [4.3 Consistency with the rest of the package](#43-consistency-with-the-rest-of-the-package)
5. [Technical design](#5-technical-design)
   - [5.1 The contract](#51-the-contract)
   - [5.2 Storage](#52-storage)
   - [5.3 Changes to loo](#53-changes-to-loo)
   - [5.4 Changes to brms](#54-changes-to-brms)
   - [5.5 Dependencies](#55-dependencies)
6. [Compatibility and lifecycle](#6-compatibility-and-lifecycle)
7. [Future plans](#7-future-plans)
8. [Testing and validation](#8-testing-and-validation)
9. [Alternatives considered](#9-alternatives-considered)
10. [Open questions](#10-open-questions)
11. [Implementation plan](#11-implementation-plan)
12. [Risks](#12-risks)

## 1. Summary

After moment matching or reloo, `loo_pred_measure()` gives correct `elpd`
only. All other measures use weights and draws that do not match. This design
adds a per-fold contract: each post-hoc method stores its weights and the
matching `ylp`, `ypred`, and `mupred` columns in the `loo` object.
`loo_pred_measure()` then replaces the matching columns of the user's matrices
with the stored columns.
It needs no knowledge of the method.

## 2. Background and motivation

A post-hoc method corrects PSIS-LOO for the observations with a high Pareto
$\hat k$. Two such methods exist today:

| Method | Package | What it changes for observation $i$ | What it discards |
| :-- | :-- | :-- | :-- |
| moment matching | loo, `loo_moment_match.default()` | `pointwise[i, ]`, `pareto_k[i]`, `psis_object$log_weights[, i]` | the transformed draws |
| reloo | brms, `reloo.brmsfit()` | `pointwise[i, ]`, sets `pareto_k[i] <- 0` | the refit `fit_j`; `psis_object` is not changed |

`loo_pred_measure()` takes `ypred`, `mupred`, and `ylp` from the original
draws. It takes the weights from `loo$psis_object`
(`R/pred_measure-compute.R`, `do_pred_measure()`). Thus:

- **Moment matching:** column $i$ of the weights refers to the transformed
  draws. Column $i$ of `ypred` refers to the original draws. They do not match.
- **Reloo:** column $i$ of the weights is still the old PSIS weight. It has a
  high $\hat k$, but the diagnostics report $\hat k = 0$.

Only `elpd` is correct, because `.elpd_pointwise()` reads it from
`loo$pointwise`. The vignette tells users to request `elpd` only.

## 3. Goals and non-goals

### Goals

- Users get correct values for all measures after any post-hoc method in §7. The first implementation covers moment matching and reloo.
- Users call `loo_pred_measure()` the same way as for plain PSIS. They never
  handle the transformed or refit draws.
- brms can write the correction through one exported loo function.
- The storage format accepts the future methods in §7 without a new format.

### Non-goals

- **`loo_subsample()`.** It does not change the draws or the weights. It
  estimates the sum from a subset of observations.
- **LOGO and K-fold in the first implementation.** The contract uses folds,
  but the implementation covers LOO only.
- **Correct values for old `loo` objects.** An object from brms 2.23.2
  `reloo()` has no stored draws. loo can only warn for it (alternative D).
- **Method-specific diagnostics in `pred_measure`.** `loo$posthoc$method`
  records the method. `pred_measure` does not use it yet.
- **`group_ids` as a fold index.** `group_ids` in `do_pred_measure()` groups
  the summaries. It does not set the held-out unit.

## 4. User-facing design (API)

### 4.1 Usage scenarios

The user asks for the predictions when they run the method. They pass the
original `posterior_predict()` and `posterior_epred()` matrices, as today.

```r
# moment matching
fit_mm <- update(fit, save_pars = brms::save_pars(all = TRUE))
loo_mm <- loo(fit_mm, moment_match = TRUE, save_psis = TRUE, save_pred = TRUE)

loo_pred_measure(
  loo = loo_mm,
  y = y,
  ypred = brms::posterior_predict(fit_mm),
  mupred = brms::posterior_epred(fit_mm),
  measure = c("elpd", "crps", "rmse")
)

# reloo: the same call
loo_re <- loo(fit, reloo = TRUE, save_psis = TRUE, save_pred = TRUE)
loo_pred_measure(loo = loo_re, y = y, ypred = ..., mupred = ..., measure = ...)
```

### 4.2 Proposed interface

**New argument `save_pred = FALSE`** in `loo_moment_match.default()`,
`loo.brmsfit()`, `reloo.brmsfit()`, and `add_criterion()`.

**New exported setter.** brms must write `loo$posthoc` from `reloo()`. brms
cannot call loo internals, so loo exports one setter. `loo_moment_match()`
uses it too.

```r
add_posthoc_draws(loo, folds, method, log_weights, ylp,
                  ypred = NULL, mupred = NULL)
```

| Argument | Type / default | Meaning |
| :-- | :-- | :-- |
| `loo` | `psis_loo` object | the object to update |
| `folds` | integer vector | the changed folds |
| `method` | character, one per fold | the post-hoc method |
| `log_weights` | $S \times$ `length(folds)` matrix | new log weights |
| `ylp` | $S \times n$ matrix | `log_lik` at the new draws; $n$ observations in `folds` |
| `ypred` | $S \times n$ matrix or `NULL` | predictions at the new draws |
| `mupred` | $S \times n$ matrix or `NULL` | expected predictions at the new draws |

Return value: the updated `loo` object.

Errors: `add_posthoc_draws()` stops if `nrow()` of a matrix is not $S$
(decision 1). It also stops if `ncol()` does not match `folds`.

**Warning cases in `loo_pred_measure()`:**

| Situation | Behaviour |
| :-- | :-- |
| `loo$posthoc` present, `ypred`/`mupred` stored | replace the columns; all measures are correct |
| `loo$posthoc` present, no stored `ypred`/`mupred`, only `elpd`-type measures requested | no message |
| `loo$posthoc` present, no stored `ypred`/`mupred`, other measures requested | warn: "Run the method again with `save_pred = TRUE`." |
| `ylp` + `psis_object` without `loo` | `psis_object` has no `posthoc`; warn if its weights differ from PSIS of `ylp` (decision 3) |

### 4.3 Consistency with the rest of the package

- `folds` matches the element name in `kfold` objects.
- `post_pred_i_upars` has the same form as the existing `log_lik_i_upars`.
- `save_pred` follows the pattern of `save_psis`.

## 5. Technical design

### 5.1 The contract

The contract does not depend on the method. It works per fold $g$, the unit
that is held out. A new element `loo$folds`, a length-$N$ vector, gives the
fold $g(i)$ of observation $i$. LOO is the case `folds = 1:N`. LOGO has one
fold per group.

For each fold $g$, the leave-out predictive distribution is one pair:

> the draws behind fold $g$, and the log weights for fold $g$.

Column $i$ of the user's matrices uses the pair of fold $g(i)$.

| Method | Draws behind fold $g$ | Log weights for fold $g$ |
| :-- | :-- | :-- |
| PSIS | full-data posterior | PSIS weights |
| moment matching | transformed draws; with `split = TRUE`, the first $S/2$ rows are $T(\theta)$ | `lwi` |
| reloo | refit on $y_{-g}$ | uniform, $-\log S$ |

**Invariant.** Column $i$ of `ylp`, `ypred`, and `mupred`, and column $g(i)$ of
`psis_object$log_weights` all refer to the same $S$ draws.

A post-hoc method must do three things to keep the invariant:

1. Write the new weights into `psis_object$log_weights[, g]`.
2. Store the matching columns of `ylp`, `ypred`, and `mupred` in the `loo`
   object.
3. Write its elpd into `pointwise[g, "elpd_loo"]`. `pointwise` has one row for
   each fold.

Rule 3 lets a method use its own elpd estimator. `pred_measure` reads `elpd`
from `loo$pointwise` (`.elpd_pointwise()`). Thus `pred_measure` needs no
change.

### 5.2 Storage

The storage has two layers.

**Base draws.** A new element `loo$draws` names the draw set that the user's
`ylp`, `ypred`, and `mupred` come from:

- `"posterior"` (default): the full-data posterior. PSIS uses it.
- `"mixture"`: the (balanced) MixIS draws. The user computes the matrices from
  the mixture fit. Its weights $\log v_{g,s}$ stay in `psis_object`.

A method that changes all folds changes the base draws. It does not use
`loo$posthoc`. For MixIS, `pred_measure` can check the user's `ylp` against the
weights: $\log v_{g,s} = \log u_s + \log\alpha_g - \text{ylp}_{g,s}$.

**Overrides.** A new element `loo$posthoc` holds only the folds that a
post-hoc method changed. It is `NULL` when no method ran:

```r
loo$posthoc <- list(
  folds  = c(12L, 57L),                # changed folds
  method = c("moment_match", "reloo"), # one entry per changed fold
  ylp    = matrix(, S, 2),             # log_lik at the new draws; one column
                                       # per observation in `folds`
  ypred  = matrix(, S, 2),             # NULL if save_pred = FALSE
  mupred = matrix(, S, 2)              # NULL if save_pred = FALSE
)
```

For each matrix, the store holds $S$ values for each observation in a changed
fold. It holds no parameter draws. Thus the size does not depend on the number
of parameters.

**Method record.** `attr(loo, "posthoc")` holds the name of each method that
changed the object, once. An example is `c("moment_match", "reloo")`. Users
can read it to see which methods the object used. `add_posthoc_draws()`
appends `method` to the attribute. `loo$posthoc$method` keeps the method for
each fold.

### 5.3 Changes to loo

| Area | Change | New / modified | Exported? |
| :-- | :-- | :-- | :-- |
| `add_posthoc_draws()` | setter, sketch below | new | yes |
| `R/loo_moment_matching.R` | `save_pred`, `post_pred_i_upars`; keep final draws | modified | `loo_moment_match()` yes |
| `R/split_moment_matching.R` | return `upars_trans_half` | modified | no |
| `R/pred_measure-compute.R`, `do_pred_measure()` | replace columns from `loo$posthoc` | modified | no |
| `R/pred_measure.R` | `@details`, `@param loo` | modified | yes |
| `NAMESPACE`, `_pkgdown.yml` | add `add_posthoc_draws` | modified | — |
| `pred-measure-workflow.Rmd` | replace the "`elpd` only" callout | modified | — |
| `NEWS.md` | new feature; note on old objects | modified | — |

**`add_posthoc_draws()`:**

```r
#' @export
add_posthoc_draws <- function(loo, folds, method, log_weights, ylp,
                              ypred = NULL, mupred = NULL) {
  # checks: nrow() of every matrix equals S, else stop (decision 1);
  # ncol(log_weights) equals length(folds); ncol() of the other matrices
  # equals the number of observations in `folds`
  loo$psis_object$log_weights[, folds] <- log_weights
  attr(loo$psis_object, "norm_const_log") <-
    matrixStats::colLogSumExps(loo$psis_object$log_weights)
  # append to loo$posthoc; a fold that is already there is replaced
  ...
  loo
}
```

**`loo_moment_match.default()`:**

- New arguments `save_pred = FALSE` and `post_pred_i_upars = NULL`.
  `post_pred_i_upars(x, upars, i, ...)` returns
  `list(ypred = <S vector>, mupred = <S vector>)`.
- If `save_pred = TRUE`, require `post_pred_i_upars`.
- `loo_moment_match_i()` keeps the final draws:
  - without split: `uparsi`;
  - with split: `upars_trans_half`. `loo_moment_match_split()` must return it.
    Today it returns only `lwi`, `lwfi`, `log_liki`, and `r_eff_i`.
- `loo_moment_match_i()` calls `post_pred_i_upars()` on the final draws. It
  then drops the draws and returns `log_liki`, `ypred_i`, and `mupred_i`.
- The update loop calls `add_posthoc_draws(loo, folds = I,
  method = "moment_match", ...)`. It replaces the direct write to
  `loo$psis_object$log_weights[, i]`.

**`do_pred_measure()`:** when `source == "loo"`, run this code before
`do_pred_measure()` reads `log_weights`:

```r
ph <- loo$posthoc %||% predperf$posthoc
if (!is.null(ph)) {
  cols <- which(loo$folds %in% ph$folds)  # equals ph$folds for LOO
  if (!is.null(ylp))    ylp[, cols]    <- ph$ylp
  if (!is.null(ypred))  ypred[, cols]  <- ph$ypred
  if (!is.null(mupred)) mupred[, cols] <- ph$mupred
  # warn if ypred/mupred were given but ph$ypred/ph$mupred are NULL
}
```

- Keep `posthoc` in the result, so that `pred_measure(predperf = ...)` can add
  measures later.
- With the replaced `ylp` columns, the `ylp` + `psis_object` pattern also
  gives the correct `elpd`.

### 5.4 Changes to brms

The code below comes from brms 2.23.2.

**Moment matching.** Add a helper next to `.log_lik_i_upars`. It uses the
existing `.update_pars()`, which returns a `brmsfit` with the transformed
draws.

```r
.post_pred_i_upars <- function(x, upars, i, newdata, resp = NULL, ...) {
  x <- update_misc_env(x, only_windows = TRUE)
  x <- .update_pars(x, upars = upars, ...)
  nd <- newdata[i, , drop = FALSE]
  list(
    ypred  = as.vector(posterior_predict(x, newdata = nd, resp = resp)),
    mupred = as.vector(posterior_epred(x, newdata = nd, resp = resp))
  )
}
```

`loo_moment_match.brmsfit()` passes `post_pred_i_upars = .post_pred_i_upars`
and `save_pred` to `loo::loo_moment_match.default()`.

**Reloo.** In `reloo.brmsfit()`, the inner function `.reloo(j)` returns only
`log_lik()` of `fit_j`. Change it to return a list:

```r
list(
  ylp    = do_call(log_lik, ll_args),
  ypred  = if (save_pred) posterior_predict(fit_j, newdata = mf[omitted, , drop = FALSE]),
  mupred = if (save_pred) posterior_epred(fit_j, newdata = mf[omitted, , drop = FALSE])
)
```

After the loop, call:

```r
loo <- loo::add_posthoc_draws(
  loo, folds = obs, method = "reloo",
  log_weights = matrix(-log(S), S, length(obs)),
  ylp = ylp, ypred = ypred, mupred = mupred
)
```

`loo$pointwise` and `pareto_k` keep their present updates.

**Pass `save_pred` through.** `loo.brmsfit()` and `add_criterion()` pass
`save_pred` to `loo_moment_match()` and `reloo()`.

### 5.5 Dependencies

No new dependencies. `matrixStats` is already used for `colLogSumExps`.

## 6. Compatibility and lifecycle

- `save_pred = FALSE` is the default. Existing calls give the same result.
- New `loo` elements (`folds`, `draws`, `posthoc`) are additive. Old objects
  lack them. `do_pred_measure()` must treat a missing `folds` as `1:N`.
- Old reloo objects have no stored draws. `loo_pred_measure()` warns for
  them (alternative D).
- loo must release `add_posthoc_draws()` before brms can call it. brms needs a
  minimum loo version.
- TODO: mark `add_posthoc_draws()` as experimental?

## 7. Future plans

The CV diagnostics workflow as worked out by Aki is consistent with the design. (Details will follow)

## 8. Testing and validation
First thoughts...
- `add_posthoc_draws()`: errors on a wrong $S$ and on a wrong `ncol()`.
  A second call on the same fold replaces the first.
- Moment matching with `save_pred = TRUE`, with and without `split`: the
  replaced columns match `post_pred_i_upars()` on the final draws.
- `do_pred_measure()`: with `posthoc`, all measures differ from the
  result without replacement only in the changed columns.
- Reloo: uniform weights give the plain mean of the refit draws.
- Snapshot tests for the three warning cases in §4.2.
- Old objects without `folds` or `posthoc` give the same result as today.

## 9. Alternatives considered

### A: Do nothing; document "`elpd` only"

The current state. No new API. Users get wrong values for other measures
without a warning, unless they read the vignette.

### B: Store the transformed or refit parameter draws

`pred_measure` could recompute the predictions. The size grows with the number
of parameters. loo cannot call `posterior_predict()` on brms objects.

Status: rejected.

### C: Full store for all folds

One layer only. For MixIS it copies the $S \times N$ matrices. With
$S = 4000$ and $N = 10\,000$, each matrix is 320 MB. The workflow also
combines the layers: MixIS for all folds, then brute-force for the folds that
fail.

Status: rejected. The two layers in Section 5.2 replace it.

### D: Warn only

`loo_pred_measure()` detects a post-hoc method. It warns if the user requests
a measure other than `elpd`, `mlpd`, or `ic`. The cost is small. But the user
cannot get correct values.

Detection:

- **Moment matching:** `loo_moment_match()` sets `attr(loo, "posthoc")` to
  `"moment_match"`. It sets the attribute only if it changes an observation.
- **Reloo:** `reloo()` sets `loo$diagnostics$pareto_k` to 0 for the refit
  observations. It does not change `loo$psis_object`. Thus the two sets of
  $\hat k$ values differ. This check also finds objects from brms 2.23.2.

*Limit:* the warning shows only when the user passes `loo` directly. A later
`pred_measure(predperf = ...)` call does not warn.

*Status:* implemented as step 0 in Section 11. Steps 1-5 replace the warning with correct values.

## 10. Open questions

| # | Question | Options | Needed by |
| :-- | :-- | :-- | :-- |
| 1 | `posterior_predict()` with one row can fail, or give wrong values, in models with `ar()`, `car()`, or `gp()`. | brms problem; add a brms test for these terms. Not part of the contract. | brms PR |
| 2 | Which diagnostics does `pred_measure` report per method ($\hat k$, MixIS trust diagnostic, MCSE)? | record `method` now; use it later | MixIS |
| 3 | Lifecycle stage of `add_posthoc_draws()`. | experimental | loo PR |

## 11. Implementation plan

Steps:

0. loo: warn only (alternative D). `loo_moment_match()` sets `attr(loo, "posthoc")`. `.warn_posthoc()` warns for measures other than `elpd`, `mlpd`, or `ic`. Status: done.
1. loo: `add_posthoc_draws()`, `loo$folds`, and the column replacement in
   `do_pred_measure()`, with tests.
2. loo: `save_pred` in `loo_moment_match.default()` and the split variant.
3. loo: documentation, vignette, `NEWS.md`.
4. brms: `.post_pred_i_upars()` and `save_pred` pass-through for moment
   matching.
5. brms: `reloo.brmsfit()` writes `posthoc`.

Steps 1-3 work without brms. Steps 4-5 need a loo release.

## 12. Risks

- **brms depends on a loo release.** Mitigation: keep the setter small and
  stable.
- **Old reloo objects give silent wrong values.** Mitigation: `NEWS.md` and documentation.
