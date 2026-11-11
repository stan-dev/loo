# Holdout predictive performance measures

Score predictive performance on **genuinely new (holdout) data** that
was not used to fit the model. This mirrors the cross-validation goal of
assessing how well a model predicts unseen observations, but with an
explicit train/test split rather than LOO or k-fold reweighting.

Supply `ylp_test` from log predictive densities evaluated on the holdout
set (e.g. `brms::log_lik(fit, newdata = test_data)`). This is required
for `elpd_test`, `mlpd_test` and `ic_test`. Optional distributional and
point-prediction measures use observed and predicted values on the test
set only. Pass training `ylp` only when an additional measure needs log
predictive densities from the training fit.

## Usage

``` r
test_pred_measure(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_test = NULL,
  measures = NULL,
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

- ylp_test:

  Matrix of pointwise log predictive densities for holdout observations
  (`S` × `n_test`), or a 3-D array (`I` × `C` × `n_test`), typically
  from `brms::log_lik(fit, newdata = test_data)`. Used with `ylp` (from
  the training fit) in `test_pred_measure()` to score genuinely new
  data.

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

- group_ids:

  Optional vector of group identifiers for grouped summaries (reserved;
  not yet implemented).

- control:

  Named list of per-measure settings. Each name must match an element of
  `measures`; the value is a list of arguments passed to that measure's
  summary function (e.g. `list(new_measure = list(add_arg = 10))`).

## Value

An object of class `"test_pred_measure"` and `"pred_measure"` with
`estimates` and `pointwise`. Measure names carry a `_test` suffix (e.g.
`elpd_test`, `rps_test`). Attribute `dims` reflects the test-set size
(from `ylp_test`), not the training data.

## Details

`elpd_test` is computed from `ylp_test` on the holdout observations
only.

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
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md),
[`kfold_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/kfold_pred_measure.md),
[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md),
[supported_measures_list](https://mc-stan.org/loo/pr/363/reference/supported_measures_list.md),
[Predictive model performance: Predictive schemes with
`pred_measure()`](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html)

## Examples

``` r
# \donttest{
if (requireNamespace("brms", quietly = TRUE)) {
  data <- lme4::sleepstudy
  train <- data[1:150, ]
  test <- data[151:nrow(data), ]
  fit <- brms::brm(
    Reaction ~ Days, data = train,
    refresh = 0, chains = 2, iter = 1000
  )
  test_pred_measure(
    y = test$Reaction,
    ypred = brms::posterior_predict(fit, newdata = test),
    mupred = brms::posterior_epred(fit, newdata = test),
    ylp = brms::log_lik(fit),
    ylp_test = brms::log_lik(fit, newdata = test),
    measures = c("rmse", "r2")
  )
}
#> Compiling Stan program...
#> Start sampling
#> 
#> Computed from 1000 posterior draws and 30 observations.
#> Data source: test
#> 
#>           Estimate    SE
#> rmse_test     26.3   2.2
#> r2_test      0.650 0.064
# }
```
