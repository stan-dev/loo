# Add predictive performance measures to an existing result

Extend a `"pred_measure"` object with additional measures **without
recomputing** what is already stored. Use this for interactive
exploration or when you first compute `elpd` and later add
distributional or point-prediction metrics.

Pass the existing object as `predperf` and supply any inputs newly
required by the requested measures (see the input table in
[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md)).
The evaluation mode (`"insample"`, `"loo"`, `"kfold"`, or `"test"`) is
taken from `predperf`; LOO paths reuse stored PSIS weights when
available.

## Usage

``` r
pred_measure(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measures = NULL,
  predperf,
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

- predperf:

  An existing predictive measure object (class `"pred_measure"`) to
  update. When supplied, stored rows and (for LOO) PSIS weights are
  reused instead of recomputed.

- group_ids:

  Optional vector of group identifiers for grouped summaries (reserved;
  not yet implemented).

- psis_object:

  A `psis` object with LOO importance weights. An alternative to passing
  a full `loo` object; must be supplied together with `ylp` when
  computing `elpd`.

- save_psis:

  Logical. If `TRUE`, store the `psis` object in the result so that
  additional measures can be added later with `pred_measure()` without
  recomputing PSIS weights.

- control:

  Named list of per-measure settings. Each name must match an element of
  `measures`; the value is a list of arguments passed to that measure's
  summary function (e.g. `list(new_measure = list(add_arg = 10))`).

## Value

An updated object of the same class as `predperf`, with new rows in
`estimates` and columns in `pointwise` for each requested measure. Rows
already in `predperf` are not recomputed. If `mlpd` or `ic` needs `elpd`
and `predperf` has no `elpd` column, `elpd` is recomputed from `ylp` for
insample and LOO results. For k-fold and test results this is an error.

## Details

**Typical workflow:**


    result <- loo_pred_measure(loo = loo_fit, save_psis = TRUE)
    pred_measure(
      y = y,
      mupred = mupred,
      predperf = result,
      measures = c("rmse", "r2")
    )

When extending a LOO result, ensure the initial call used
`save_psis = TRUE` (or that `predperf` already contains a `psis_object`)
so LOO weights are available for additional measures.

## See also

[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md),
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md),
[`kfold_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/kfold_pred_measure.md),
[`test_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/test_pred_measure.md),
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
  result <- insample_pred_measure(
    ylp = brms::log_lik(fit),
    y = fit$data$Reaction,
    ypred = brms::posterior_predict(fit)
  )
  pred_measure(
    y = fit$data$Reaction,
    mupred = brms::posterior_epred(fit),
    predperf = result,
    measures = c("rmse", "r2")
  )
}
#> Compiling Stan program...
#> Start sampling
#> 
#> Computed from 1000 posterior draws and 180 observations.
#> Data source: in-sample
#> 
#>      Estimate     SE
#> elpd   -950.1   10.1
#> rmse     47.4    2.7
#> r2      0.286  0.057
# }
```
