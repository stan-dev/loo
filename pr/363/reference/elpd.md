# Generic (expected) log-predictive density (deprecated)

As of loo 3.0.0, `elpd()` is **deprecated**. Please use
[`measure_elpd()`](https://mc-stan.org/loo/pr/363/reference/measure_elpd.md)
instead. For full predictive performance workflows, see
[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md)
and
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md).
See
[`vignette("migration-guide", package = "loo")`](https://mc-stan.org/loo/pr/363/articles/migration-guide.md)
for a full mapping table.

## Usage

``` r
elpd(x, ...)

# S3 method for class 'array'
elpd(x, ...)

# S3 method for class 'matrix'
elpd(x, ...)
```

## Arguments

- x:

  A log-likelihood array or matrix. The **Methods (by class)** section,
  below, has detailed descriptions of how to specify the inputs for each
  method.

- ...:

  Currently ignored.

## Details

The return type differs: `elpd()` returns class `"elpd_generic"` with
`elpd` and `ic` in `pointwise`;
[`measure_elpd()`](https://mc-stan.org/loo/pr/363/reference/measure_elpd.md)
returns class `"measure"`.

The `elpd()` methods for arrays and matrices can compute the expected
log pointwise predictive density for a new dataset or the log pointwise
predictive density of the observed data (an overestimate of the elpd).
The `elpd()` function is an S3 generic and methods are provided for 3-D
pointwise log-likelihood arrays and matrices.

## Methods (by class)

- `elpd(array)`: An \\I\\ by \\C\\ by \\N\\ array, where \\I\\ is the
  number of MCMC iterations per chain, \\C\\ is the number of chains,
  and \\N\\ is the number of data points.

- `elpd(matrix)`: An \\S\\ by \\N\\ matrix, where \\S\\ is the size of
  the posterior sample (with all chains merged) and \\N\\ is the number
  of data points.

## See also

[`measure_elpd()`](https://mc-stan.org/loo/pr/363/reference/measure_elpd.md),
[`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md),
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md),
and the vignette *Holdout validation and K-fold cross-validation*.

## Examples

``` r
# \dontrun{
# Deprecated:
LLarr <- example_loglik_array()
elpd(LLarr)
#> Warning: 'elpd.array' is deprecated.
#> Use 'measure_elpd' instead.
#> See help("Deprecated")
#> 
#> Computed from 1000 posterior draws and 32 log-likelihood terms using the generic elpd function.
#> 
#>      Estimate  SE
#> elpd    -80.3 3.2
#> ic      160.5 6.5
# ->
measure_elpd(LLarr)
#> 
#> Computed from 1000 draws by 32 observations.
#> 
#>      Estimate    SE
#> elpd    -80.3   3.2
# }
```
