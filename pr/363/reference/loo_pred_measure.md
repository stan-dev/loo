# PSIS-LOO predictive performance measures

Estimate out-of-sample predictive performance with **PSIS-LOO**
(Pareto-smoothed importance sampling leave-one-out cross-validation).
PSIS-LOO approximates exact LOO-CV without refitting the model once per
observation: each held-out point is scored by reweighting the full-data
posterior draws.

With the default `measures = NULL`, the result holds `elpd_loo`, the LOO
estimate of expected log pointwise predictive density (ELPD), and
`p_loo` (effective number of parameters, the difference between
in-sample and LOO log predictive density). `p_loo` is reported together
with `elpd_loo`. See
[`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md) and the
[Cross-validation FAQ](https://users.aalto.fi/~ave/CV-FAQ.html) for
interpretation.

## Usage

``` r
loo_pred_measure(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measures = NULL,
  loo = NULL,
  group_ids = NULL,
  psis_object = NULL,
  save_psis = FALSE,
  control = list()
)
```

## Arguments

- y:

  Vector of observed values (`n`). Required for every measure except
  `elpd`, `mlpd`, and `ic`.

- ypred:

  Matrix of posterior predictive draws (`S` draws × `n` observations),
  typically from
  [`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html).
  Required for `rps`, `srps`, and `brier`.

- mupred:

  Matrix of posterior expected values (`S` × `n`), typically from
  [`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html).
  Required for `mae`, `mse`, `rmse`, `r2`, `acc`, and `bacc`.

- ylp:

  Matrix of pointwise log predictive densities or probabilities (`S` ×
  `n`), typically from
  [`brms::log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  or a 3-D array (`I` × `C` × `n`: iterations, chains, observations).
  Required for density-based summaries (`elpd`, `mlpd`, `ic`).

- measures:

  Measures to compute. `NULL` (default) gives `elpd` (with `p_loo` /
  `p_kfold` for LOO and k-fold). Can be:

  - A **character vector** of built-in names; see
    [supported_measures_list](https://mc-stan.org/loo/pr/363/reference/supported_measures_list.md).

  - A **function** built with
    [`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md)
    for one custom measure.

  - A **list** mixing character scalars (built-in names) and named
    functions (custom measures), e.g. `list("rps", my_metric = my_fun)`.

  Custom functions are called with any of `y`, `ypred`, `mupred`, `ylp`,
  and `log_weights` that appear in their formals, plus arguments from
  `control`. They must return a list with `estimate`, `se`, and
  `pointwise`. A length-2 `estimates` (estimate and SE) can replace
  `estimate` and `se`.

  A custom measure declares whether it is a loss (lower is better) or a
  utility (higher is better) with `loss` in
  [`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md).
  Without it a custom measure is taken to be a utility.
  [`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
  uses the declaration to put all measures on a common utility scale and
  to rank models, so an undeclared loss is compared and ranked in the
  wrong direction.

  A custom measure can declare how the standard error of a difference
  between two models is computed, with `se_diff_fun` in
  [`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md).
  It accepts a function `function(ref, cmp) ...` or the shorthands
  `"sum"` and `"mean"` for the paired pointwise formulas. Without a
  declaration,
  [`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
  reports the difference with an `NA` standard error. A function
  receives one list per model with elements `estimate`, `se`,
  `pointwise`, and `extra`, always on the measure's natural scale, and
  must return the standard error of the difference as a numeric scalar.

  `extra` is for anything the standard error needs that the pointwise
  values do not carry. Return it as an additional list element `extra`
  from the measure function and it is stored alongside the estimates and
  passed on to its `measure_se_diff` function; the built-in `r2` uses it
  for the baseline `(y_i - mean(y))^2`, which cannot be recovered once
  `y` is out of scope.

- loo:

  A [`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md) result,
  computed with `save_psis = TRUE` so that PSIS weights are available
  for additional measures. See `loo_pred_measure()`.

- group_ids:

  Optional vector of group identifiers for grouped summaries (reserved;
  not yet implemented).

- psis_object:

  A `psis` object with LOO importance weights. An alternative to passing
  a full `loo` object; must be supplied together with `ylp` when
  computing `elpd`.

- save_psis:

  Logical. If `TRUE`, store the `psis` object in the result so that
  additional measures can be added later with
  [`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
  without recomputing PSIS weights.

- control:

  Named list of per-measure settings. Each name must match an element of
  `measures`; the value is a list of arguments passed to that measure's
  summary function (e.g. `list(new_measure = list(add_arg = 10))`).

## Value

An object of class `"loo_pred_measure"`, `"pred_measure"`, and `"loo"`.
When a `loo` object is supplied, additional classes from that object are
inherited (e.g. `"psis_loo"`). In addition to `estimates` and
`pointwise`, the list may contain:

- `diagnostics`:

  PSIS diagnostics, including Pareto \\\hat{k}\\ in
  `diagnostics$pareto_k`. A value above the sample-size dependent
  threshold \\\min(1 - 1 / \log\_{10}(S), 0.7)\\, where \\S\\ is the
  number of draws, suggests an unreliable LOO estimate for that
  observation. The [`print()`](https://rdrr.io/r/base/print.html) method
  reports the threshold it used.

- `log_weights`:

  Normalized log importance weights used for LOO scoring.

- `psis_object`:

  Stored when `save_psis = TRUE`; needed to add further measures with
  [`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
  without recomputing weights.

Measure names carry a `_loo` suffix (e.g. `elpd_loo`, `rps_loo`).

## Details

**Three input patterns:**

- Precomputed `loo` object:

  `loo_pred_measure(loo = loo_fit, ...)`. Run
  [`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md) with
  `save_psis = TRUE`.

- `ylp` + `psis_object`:

  Pass both when you have already computed PSIS weights separately.

- `ylp` only:

  PSIS weights are computed internally from `ylp`.

If you corrected the `loo` object after PSIS, pass it as `loo`. For
example,
[`loo_moment_match()`](https://mc-stan.org/loo/pr/363/reference/loo_moment_match.md)
or
[`brms::reloo()`](https://paulbuerkner.com/brms/reference/reloo.brmsfit.html).

For distributional and point-prediction measures (`rps`, `r2`, etc.),
supply `y`, `ypred`, and/or `mupred` as for
[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md).
When adding measures incrementally, call
[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
with `predperf` set to an existing result; use `save_psis = TRUE` on the
initial call so weights are stored.

## `measure_info` attribute

All `*_pred_measure()` and
[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
results include attribute `measure_info`: a named list, keyed by bare
measure name, recording what
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
needs to know about each measure — `loss` (whether lower values are
better), `diff_method` (how the standard error of a difference is
obtained), and, where applicable, `se_diff_fun` and `extra`. When
measures are added incrementally with
[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md),
the attribute is extended for the newly computed measures.

Built-in measures take `loss`, `diff_method`, and `se_diff_fun` from the
package measure registry. Custom measures always get
`diff_method = "custom"` and take `loss` and `se_diff_fun` from
[`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md).
Without `loss = TRUE` a custom measure is treated as a utility, so an
undeclared loss is compared and ranked in the wrong direction. See
[loo-glossary](https://mc-stan.org/loo/pr/363/reference/loo-glossary.md).

## See also

[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md),
[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md),
[`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md),
[supported_measures_list](https://mc-stan.org/loo/pr/363/reference/supported_measures_list.md),
[Predictive model performance: Predictive schemes with
`pred_measure()`](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html)

## Examples

``` r
# \donttest{
if (requireNamespace("brms", quietly = TRUE)) {
  fit <- brms::brm(
    Reaction ~ Days, data = lme4::sleepstudy,
    refresh = 0, chains = 2, iter = 1000
  )
  loo_fit <- loo::loo(fit, save_psis = TRUE)
  loo_pred_measure(
    loo = loo_fit,
    y = fit$data$Reaction,
    ypred = brms::posterior_predict(fit),
    measures = c("rmse", "r2")
  )
}
#> Compiling Stan program...
#> Start sampling
#> Error in .validate_numeric_matrix(mupred, arg = "mupred", ncol = length(y)): `mupred` must be a numeric matrix or array, not NULL.
# }
```
