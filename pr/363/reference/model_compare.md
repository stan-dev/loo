# Model comparison

Compare fitted models on
[ELPD](https://mc-stan.org/loo/pr/363/reference/loo-glossary.md) or, for
[`pred_measure`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
results, on several predictive performance measures at once.

`model_compare()` accepts two families of input:

- **Classic results** — `"loo"`, `"waic"`, and `"kfold"` objects,
  compared on ELPD alone.

- **Predictive measure results** — objects from
  [`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md),
  [`kfold_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/kfold_pred_measure.md),
  [`test_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/test_pred_measure.md),
  or
  [`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md),
  compared on every measure the models share.

All models in one call must be evaluated the same way. Differences
between, say, a LOO and a k-fold result would contrast held-out schemes
rather than models, so mixed inputs are an error.

## Usage

``` r
# S3 method for class 'compare.loo'
print(x, ..., digits = NULL, p_worse = TRUE, simplify = TRUE, measures = NULL)

model_compare(x, ...)

# Default S3 method
model_compare(x, ...)
```

## Arguments

- x:

  An object of class `"loo"` or `"pred_measure"`, or a list of such
  objects. List names are used as the model names in the output. See
  **Examples**.

- ...:

  Additional objects of class `"loo"` or `"pred_measure"`, if not passed
  in as a single list. Naming every model here, as in
  `model_compare(A = m1, B = m2)`, names the models in the output,
  exactly as the list form does.

- digits:

  For the print method only, the number of decimal places to print. The
  default `NULL` gives each measure its own: 1 for `elpd` and `ic`, and
  3 for `mlpd` and for the bounded measures (`r2`, `acc`, `bacc`,
  `brier`). A measure on the scale of the data (`mae`, `rmse`, `mse`,
  `rps`, `srps`, and any custom measure) takes the places that show two
  significant digits of its own standard error, at most 4. Pass a single
  number for one format in every column, or a named vector such as
  `c(r2 = 2, rmse = 3)` to override single measures. `p_worse` always
  prints with 2 places.

- p_worse:

  For the print method only, should we include the normal approximation
  based probability of each model having worse performance than the
  reference model? The default is `TRUE`.

- simplify:

  For the print method only, should the output be simplified to only
  include the model names, differences, and (when `p_worse = TRUE`)
  diagnostic columns? The default is `TRUE`. Set to `FALSE` to also
  print the available estimate columns: pointwise ELPD, LOOIC/WAIC and
  their standard errors for classic comparisons. For
  [`pred_measure`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
  comparisons each printed table gains the estimate and standard error
  of its own measure, and an ELPD table also gains `p` and `se_p`. With
  the default `measures = NULL` one table is printed, and it gains the
  estimates of every compared measure. The difference columns are never
  added.

- measures:

  For `pred_measure` comparisons only, which measures to print diff
  tables for. `NULL` (default) prints only the ranking measure (the
  first compared measure, `"elpd"` when present); `"all"` prints all
  compared measures; or a character vector of measure names (e.g.
  `c("elpd", "mse")`). Each table is sorted by its own measure, best
  model first, so the same model need not lead every table.

## Value

A data frame of class `"compare.loo"` with one row per model and its own
print method.

For classic `"loo"` / `"waic"` / `"kfold"` comparisons the columns are
unchanged from previous versions: `model`, `elpd_diff`, `se_diff`,
`p_worse`, `diag_diff`, `diag_elpd`, and the estimate columns of the
input objects.

For
[`pred_measure`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
comparisons there is a `{measure}_diff` and a `{measure}_se_diff` column
for every measure shared by all models (e.g. `rmse_diff`,
`rmse_se_diff`). ELPD-family measures use `elpd_diff` and `se_diff`
instead. `p_worse` and `diag_diff` are computed for ELPD only.
`diag_elpd` holds per-model Pareto \\\hat{k}\\ diagnostics and is
present only for
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)
comparisons, the only source with Pareto \\\hat{k}\\ values.

The object also carries the following attributes:

- `compare_reference`:

  A named character vector giving, for each measure, the model its
  differences were computed against, which is that measure's own best
  model.

- `compare_measures`:

  Bare names of all measures that were compared.

- `sign_converted_measures`:

  Bare names of the loss measures whose sign was flipped onto the
  utility scale.

- `compare_source`:

  The shared evaluation source: `"loo"`, `"kfold"`, `"test"`, or
  `"insample"`.

`compare_reference` is set for every comparison; the last three are set
for
[`pred_measure`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
comparisons only.

## Details

### Differences and their standard errors

Differences are pairwise: every model is compared with one reference
model, whose own `{measure}_diff` is therefore `0`. The reference is the
best model on that measure, so `mse_diff` may use a different reference
than `elpd_diff`, and the remaining differences for a measure are all
negative. Rows are ordered by `"elpd"` when all models share it.
Otherwise, rows are ordered by the first shared measure in alphabetical
order.

The standard error of a difference is a paired estimate, which uses the
fact that the same \\N\\ data points were used for both models. It
should not be expected to equal the difference of the two models'
standard errors.

### `p_worse`, `diag_diff`, and `diag_elpd`

`p_worse` is the probability that a model has worse ELPD than the
reference model, computed with a normal approximation from `elpd_diff`
and `se_diff`. Sivula et al. (2025) give the conditions under which that
approximation is good; `diag_diff` reports the two that fail most often:

- `N < 100` (small data)

- `|elpd_diff| < 4` (models make similar predictions)

Either message means the error distribution is skewed or thick tailed,
the normal approximation is not well calibrated, and `p_worse` is likely
too large. If `|elpd_diff|` is many times `se_diff` the difference is
quite certain. Model misspecification and outliers also skew the error
distribution, and can be diagnosed with the usual predictive checks.

`diag_elpd` reports the PSIS-LOO Pareto \\\hat{k}\\ diagnostic for each
model's pointwise ELPD. An entry `K k_psis > 0.7`, where `K` counts the
high Pareto \\\hat{k}\\ values, warns of possible bias in `elpd_diff`
favoring models with many such values. Pareto \\\hat{k}\\ describes a
model's PSIS-LOO approximation rather than any one measure or pair of
models, and every LOO measure uses the same importance weights, so for
`pred_measure` comparisons
[`print()`](https://rdrr.io/r/base/print.html) reports it once per model
in a block above the difference tables instead of as a column inside one
of them. The `diag_elpd` column is still returned on the object.

### Comparing `pred_measure` objects

When all inputs are predictive measure results sharing one evaluation
source, paired differences are computed for every measure present in all
models. Measures are matched on their bare names, so the source suffix
(`_loo`, `_kfold`, `_test`, or none for in-sample) is handled
transparently. When the models were evaluated with different `measures`,
only the shared measures are compared and a warning lists the omitted
ones.

The data frame carries one row order for all measures, but each
*printed* measure table is sorted by its own difference, so the best
model on that measure always leads its table and the differences run in
decreasing order. Use `print(x, measures = "all")` to display a table
for every compared measure; see
[loo-glossary](https://mc-stan.org/loo/pr/363/reference/loo-glossary.md)
for column definitions.

### Utility scale and sign conversion

Measures differ in orientation in their raw form: ELPD and SRPS/SCRPS
are utilities (higher is better), while MSE, RPS/CRPS and the Brier
score are losses (lower is better). All `{measure}_diff` values are
reported on a common utility scale, so loss measures have their sign
flipped and a negative `{measure}_diff` always means worse performance
than the reference. Which measures are losses is recorded in the `loss`
element of each measure's entry in the `measure_info` attribute of an
`*_pred_measure()` result. The flipped measures are named in the
`sign_converted_measures` attribute.
[`print()`](https://rdrr.io/r/base/print.html) marks them with "sign
flipped" in the table header and names them below the tables.

A custom measure is treated as a utility unless it declares otherwise
with `loss = TRUE` in
[`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md).
The declaration also determines the ranking direction, so an undeclared
loss is both flipped and ranked in the wrong direction; see
[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md).

### Standard error of a measure difference

How `{measure}_se_diff` is obtained is recorded in the `diff_method`
element of the measure's entry in `measure_info`:

- `"sum"` or `"mean"`: the overall estimate is the sum (`elpd`, `ic`) or
  the mean (`mlpd`, `mae`, `mse`, `acc`, `rps`, `srps`, `brier`) of its
  pointwise contributions, so the standard error is computed from paired
  pointwise differences (the same formula as `se_diff`).

- `"measure_specific"`: the overall estimate is not a sum or mean of
  pointwise contributions (`r2`, `rmse`, `bacc`), so the measure
  supplies its own standard error of the difference.

- `"custom"`: the standard error comes from what the user passes to
  [`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md).
  `{measure}_se_diff` is `NA` when the measure declares nothing.

### Source-specific behavior

Comparisons behave the same way across sources, with three exceptions:

- **`diag_elpd`** is produced only for
  [`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)
  comparisons, since Pareto \\\hat{k}\\ diagnostics exist only for
  PSIS-LOO.

- **K-fold** comparisons warn when the models do not share the same
  number of folds, matching the behavior for plain `"kfold"` objects.

- **In-sample** comparisons warn that in-sample scores are
  optimistically biased and favor more complex models. They are
  supported for completeness, but out-of-sample sources should be
  preferred for model selection.

### Warnings for many model comparisons

If more than \\11\\ models are compared, we internally recompute the
model differences using the median model (by ELPD, or by the ranking
measure for `pred_measure` comparisons) as the baseline, and estimate
whether the differences in predictive performance are potentially due to
chance as described by McLatchie and Vehtari (2023). This flags a
warning if there is a risk of over-fitting due to the selection process.
In that case users are recommended to avoid model selection based on
LOO-CV, and instead to favor model averaging/stacking or projection
predictive inference.

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4
([journal
version](https://link.springer.com/article/10.1007/s11222-016-9696-4),
[preprint arXiv:1507.04544](https://arxiv.org/abs/1507.04544)).

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72):1-58. [PDF](https://jmlr.org/papers/v25/19-556.html)

Sivula, T, Magnusson, M., Matamoros A. A., and Vehtari, A. (2025).
Uncertainty in Bayesian leave-one-out cross-validation based model
comparison. *Bayesian Analysis*.
[doi:10.1214/25-BA1569](https://doi.org/10.1214/25-BA1569)

McLatchie, Y., and Vehtari, A. (2024). Efficient estimation and
correction of selection-induced bias with order statistics. *Statistics
and Computing*. 34(132).
[doi:10.1007/s11222-024-10442-4](https://doi.org/10.1007/s11222-024-10442-4)

## See also

- The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html)
  on the **loo** website for answers to frequently asked questions.

- The article [Model comparison: Explanation of
  `model_compare()`](https://mc-stan.org/loo/articles/articles-online-only/model-comparison.html)
  on the **loo** website, for how the differences and their standard
  errors are computed for each measure.

## Examples

``` r
# very artificial example, just for demonstration!
LL <- example_loglik_array()
loo1 <- loo(LL)     # should be worst model when compared
loo2 <- loo(LL + 1) # should be second best model when compared
loo3 <- loo(LL + 2) # should be best model when compared

comp <- model_compare(loo1, loo2, loo3)
print(comp, digits = 2)
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model3      0.00    0.00      NA                    
#>  model2    -32.00    0.00    1.00   N < 100          
#>  model1    -64.00    0.00    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.

# can use a list of objects with custom names
# the names will be used in the output
model_compare(list("apple" = loo1, "banana" = loo2, "cherry" = loo3))
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  cherry       0.0     0.0      NA                    
#>  banana     -32.0     0.0    1.00   N < 100          
#>   apple     -64.0     0.0    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.

# \dontrun{
# works for waic (and kfold) too
model_compare(waic(LL), waic(LL - 10))
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model1       0.0     0.0      NA                    
#>  model2    -320.0     0.0    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.

# compare multiple predictive measures from loo_pred_measure()
if (requireNamespace("brms", quietly = TRUE)) {
  fit1 <- brms::brm(
    Reaction ~ Days, data = lme4::sleepstudy,
    refresh = 0, chains = 2, iter = 1000
  )
  fit2 <- brms::brm(
    Reaction ~ poly(Days, 2), data = lme4::sleepstudy,
    refresh = 0, chains = 2, iter = 1000
  )
  pm1 <- loo_pred_measure(
    loo = loo(fit1, save_psis = TRUE),
    y = fit1$data$Reaction,
    mupred = brms::posterior_epred(fit1),
    measures = c("rmse", "r2")
  )
  pm2 <- loo_pred_measure(
    loo = loo(fit2, save_psis = TRUE),
    y = fit2$data$Reaction,
    mupred = brms::posterior_epred(fit2),
    measures = c("rmse", "r2")
  )
  comp <- model_compare(pm1, pm2)
  print(comp)
  print(comp, measures = "all")

  # the same works for k-fold CV
  folds <- kfold_split_random(K = 5, N = nrow(lme4::sleepstudy))
  kf1 <- brms::kfold(fit1, folds = folds, save_fits = TRUE)
  kf2 <- brms::kfold(fit2, folds = folds, save_fits = TRUE)
  kpm1 <- kfold_pred_measure(
    y = fit1$data$Reaction,
    mupred = brms::kfold_predict(kf1, method = "fitted")$yrep,
    kfold = kf1,
    measures = "rmse"
  )
  kpm2 <- kfold_pred_measure(
    y = fit2$data$Reaction,
    mupred = brms::kfold_predict(kf2, method = "fitted")$yrep,
    kfold = kf2,
    measures = "rmse"
  )
  model_compare(kpm1, kpm2)

  # mixing evaluation sources is an error
  try(model_compare(pm1, kpm2))
}
#> Compiling Stan program...
#> Start sampling
#> Compiling Stan program...
#> Start sampling
#> Each measure compared against its own best model (r2: fit1, rmse: fit1).
#>  model r2_diff r2_se_diff
#>   fit1   0.000      0.000
#>   fit2  -0.005      0.006
#> Use print(x, measures = "all") to see all measures.
#> Each measure compared against its own best model (r2: fit1, rmse: fit1).
#> 
#> -- r2 (vs fit1) --
#>  model r2_diff r2_se_diff
#>   fit1   0.000      0.000
#>   fit2  -0.005      0.006
#> 
#> -- rmse (vs fit1, sign flipped) --
#>  model rmse_diff rmse_se_diff
#>   fit1      0.00         0.00
#>   fit2     -0.16         0.20
#> 
#> All differences: 0 = best model, negative = worse.
#> Signs flipped for loss measures: rmse.
#> Setting 'K' to the number of folds (5)
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
#> Setting 'K' to the number of folds (5)
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
#> Error : All models must be evaluated on the same source, but got: loo, k-fold. Recompute all models with the same *_pred_measure() function.
# }
```
