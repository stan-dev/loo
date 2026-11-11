# Shared parameters for predictive measure wrappers

Parameter definitions shared by the user-facing entry points and the
internal engine `do_pred_measure()`.

## Usage

``` r
do_pred_measure(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_test = NULL,
  measures = NULL,
  predperf = NULL,
  loo = NULL,
  kfold = NULL,
  group_ids = NULL,
  psis_object = NULL,
  save_psis = FALSE,
  source = NULL,
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

- ylp_test:

  Matrix of pointwise log predictive densities for holdout observations
  (`S` × `n_test`), or a 3-D array (`I` × `C` × `n_test`), typically
  from `brms::log_lik(fit, newdata = test_data)`. Used with `ylp` (from
  the training fit) in
  [`test_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/test_pred_measure.md)
  to score genuinely new data.

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

- predperf:

  An existing predictive measure object (class `"pred_measure"`) to
  update. When supplied, stored rows and (for LOO) PSIS weights are
  reused instead of recomputed.

- loo:

  A [`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md) result,
  computed with `save_psis = TRUE` so that PSIS weights are available
  for additional measures. See
  [`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md).

- kfold:

  A `kfold` object from
  [`brms::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html).
  Supplies ELPD summaries and fold structure for
  [`kfold_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/kfold_pred_measure.md);
  pass `y`, `ypred`, and/or `mupred` when requesting additional
  measures.

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

- source:

  Character string indicating the evaluation mode: `"insample"`,
  `"loo"`, `"kfold"`, or `"test"`. Set automatically by the wrapper
  functions; required when calling `do_pred_measure()` directly.

- control:

  Named list of per-measure settings. Each name must match an element of
  `measures`; the value is a list of arguments passed to that measure's
  summary function (e.g. `list(new_measure = list(add_arg = 10))`).
