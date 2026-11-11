# In-sample predictive performance measures

Compute predictive performance measures on the same data used to fit the
model. This is the simplest entry point when you want density scores
(`elpd`) and optional distributional or point-prediction metrics in one
call.

In-sample `elpd` sums the expected log pointwise predictive density
(ELPD) over the training observations. Because the model has already
seen these data, in-sample scores are **optimistically biased** for
predicting future or otherwise unseen observations. For out-of-sample
performance, use
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md),
[`kfold_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/kfold_pred_measure.md),
or
[`test_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/test_pred_measure.md)
instead; see Vehtari's [Cross-validation
FAQ](https://users.aalto.fi/~ave/CV-FAQ.html) (Section 4).

## Usage

``` r
insample_pred_measure(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measures = NULL,
  group_ids = NULL,
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

- group_ids:

  Optional vector of group identifiers for grouped summaries (reserved;
  not yet implemented).

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

An object of class `"insample_pred_measure"` and `"pred_measure"`: a
list with:

- `estimates`:

  Matrix of summary estimates and standard errors (rows are measures,
  columns are `Estimate` and `SE`). The row `elpd` is present when
  `measures` is `NULL` (default) or names "elpd".

- `pointwise`:

  Matrix of observation-level contributions (one column per measure).

The attribute `source` is `"insample"`. Attribute `dims` gives posterior
draws × observations. Attribute `measure_info` records what
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
needs to know about each measure; see section below. Use
[`print()`](https://rdrr.io/r/base/print.html) for a readable summary
table.

## Details

**Input requirements by measure.** Supply only the inputs each measure
needs:

|                            |       |     |         |          |
|----------------------------|-------|-----|---------|----------|
| Measure                    | `ylp` | `y` | `ypred` | `mupred` |
| `elpd`, `mlpd`, `ic`       | ✓     |     |         |          |
| `rps`, `srps`, `brier`     |       | ✓   | ✓       |          |
| `acc`, `bacc`              |       | ✓   |         | ✓        |
| `mae`, `mse`, `rmse`, `r2` |       | ✓   |         | ✓        |

"elpd" is the default measure. `mlpd` and `ic` are derived from `elpd`,
but requesting them does not add an `elpd` row. Request `ic`, `mlpd`, or
other density scores via `measures`, or supply a custom function; see
[supported_measures_list](https://mc-stan.org/loo/pr/363/reference/supported_measures_list.md)
and the [Predictive model performance: Overview of
measures](https://mc-stan.org/loo/articles/articles-online-only/overview-measures.html)
article for definitions and orientation (higher vs lower is better).

**Custom measures.** Build a custom measure with
[`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md).
The function must return `estimate`, `se`, and `pointwise`. Only
arguments declared in the function signature among `y`, `ypred`,
`mupred`, `ylp`, and `log_weights` are supplied automatically.

Custom measures are assumed to be on a utility scale (higher is better)
in
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md).
Declare a custom loss with `loss = TRUE` in
[`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md)
so that
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
converts and ranks it in the right direction.

Declare how the standard error of a difference between two models is
computed with `se_diff_fun` in
[`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md).
It accepts a function, `"sum"`, or `"mean"`. For a measure that declares
nothing,
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
reports an `NA` standard error.

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

[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
to add measures incrementally,
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
  insample_pred_measure(
    ylp = brms::log_lik(fit),
    y = fit$data$Reaction,
    ypred = brms::posterior_predict(fit),
    mupred = brms::posterior_epred(fit),
    measures = c("rmse", "r2")
  )
}
#> Compiling Stan program...
#> Start sampling
#> 
#> Computed from 1000 posterior draws and 180 observations.
#> Data source: in-sample
#> 
#>      Estimate    SE
#> rmse     47.4   2.7
#> r2      0.286 0.057
# }
# \dontrun{
# Custom measure (same contract as built-in point-prediction metrics)
my_abs_err <- function(y, mupred, log_weights = NULL) {
  mu <- colMeans(mupred)
  pw <- abs(y - mu)
  list(
    estimate = mean(pw),
    se = sd(pw) / sqrt(length(pw)),
    pointwise = pw
  )
}
# the estimate is the mean of the pointwise values, so declare "mean"
my_abs_err <- custom_measure(my_abs_err, name = "my_abs_err", se_diff_fun = "mean")
# insample_pred_measure(y = y, mupred = mupred, ylp = ylp, measures = my_abs_err)
# }
```
