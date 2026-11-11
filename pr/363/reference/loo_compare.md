# Model comparison (deprecated)

**This function is deprecated**. Please use the new
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
function instead. See
[`vignette("migration-guide", package = "loo")`](https://mc-stan.org/loo/pr/363/articles/migration-guide.md)
for details.

`loo_compare()` remains an exported generic so that `loo_compare`
methods registered by other packages keep dispatching as before, but it
is frozen at its previous behavior: it compares `"loo"`, `"waic"`, and
`"kfold"` objects on
[ELPD](https://mc-stan.org/loo/pr/363/reference/loo-glossary.md) only.
Comparing
[`pred_measure`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
results requires
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md).

The deprecation warning is issued once per session, so a script that
calls `loo_compare()` repeatedly is not flooded with warnings.

## Usage

``` r
loo_compare(x, ...)

# Default S3 method
loo_compare(x, ...)

# S3 method for class 'psis_loo_ss_list'
loo_compare(x, ...)
```

## Arguments

- x:

  An object of class `"loo"` or a list of such objects. If a list is
  used then the list names will be used as the model names in the
  output.

- ...:

  Additional objects of class `"loo"`, if not passed in as a single
  list.

## Value

See
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md).
For the inputs `loo_compare()` still accepts, the result is identical to
what
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
returns.

## See also

[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)

## Examples

``` r
LL <- example_loglik_array()
loo1 <- loo(LL)
loo2 <- loo(LL + 1)

# deprecated; identical to model_compare(loo1, loo2)
suppressWarnings(loo_compare(loo1, loo2))
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model2       0.0     0.0      NA                    
#>  model1     -32.0     0.0    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.

# use this instead
model_compare(loo1, loo2)
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model2       0.0     0.0      NA                    
#>  model1     -32.0     0.0    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```
