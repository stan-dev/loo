# K-fold cross-validation predictive performance measures

Compute predictive performance measures under **k-fold
cross-validation**. K-fold CV holds out groups of observations, refits
(or reuses stored fits), and scores the held-out folds.

Pass a `kfold` object from
[`brms::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html)
(with `save_fits = TRUE` when you need posterior predictions on held-out
folds). `elpd_kfold` and `p_kfold` come from the `kfold` object and are
reported by default or when `measures` names "elpd"; additional measures
require the same optional inputs as
[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md).

## Usage

``` r
kfold_pred_measure(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measures = NULL,
  kfold = NULL,
  group_ids = NULL,
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

- kfold:

  A `kfold` object from
  [`brms::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html).
  Supplies ELPD summaries and fold structure for `kfold_pred_measure()`;
  pass `y`, `ypred`, and/or `mupred` when requesting additional
  measures.

- group_ids:

  Optional vector of group identifiers for grouped summaries (reserved;
  not yet implemented).

- control:

  Named list of per-measure settings. Each name must match an element of
  `measures`; the value is a list of arguments passed to that measure's
  summary function (e.g. `list(new_measure = list(add_arg = 10))`).

## Value

An object of class `"kfold_pred_measure"` and `"pred_measure"`,
inheriting attributes from the `kfold` object (`K`, `folds`,
`fold_type`, etc.). The list contains `estimates` and `pointwise`;
measure names carry a `_kfold` suffix (e.g. `elpd_kfold`, `rps_kfold`).

## Details

For distributional measures on held-out folds, obtain posterior
predictions with
[`brms::kfold_predict()`](https://paulbuerkner.com/brms/reference/kfold_predict.html)
and pass the resulting `yrep` matrices as `ypred` and/or `mupred`. See
the sleep-study workflow in [Predictive model performance: Predictive
schemes with
`pred_measure()`](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html).

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

[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md),
[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md),
[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md),
[`brms::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html),
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
  kf <- brms::kfold(fit, K = 5, save_fits = TRUE)
  ypred_kf <- brms::kfold_predict(kf, method = "predict")$yrep
  mupred_kf <- brms::kfold_predict(kf, method = "fitted")$yrep
  kfold_pred_measure(
    y = fit$data$Reaction,
    ypred = ypred_kf,
    mupred = mupred_kf,
    kfold = kf,
    measures = "rmse"
  )
}
#> Compiling Stan program...
#> Start sampling
#> Fitting model 1 out of 5
#> Start sampling
#> Fitting model 2 out of 5
#> Start sampling
#> Fitting model 3 out of 5
#> Start sampling
#> Fitting model 4 out of 5
#> Start sampling
#> Fitting model 5 out of 5
#> Start sampling
#> 
#> Computed from 1000 posterior draws and 180 observations.
#> Data source: k-fold
#> 
#>            Estimate   SE
#> rmse_kfold     48.1  2.8
# }
```
