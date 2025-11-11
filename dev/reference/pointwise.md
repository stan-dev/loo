# Convenience function for extracting pointwise estimates

Convenience function for extracting pointwise estimates

## Usage

``` r
pointwise(x, estimate, ...)

# S3 method for class 'loo'
pointwise(x, estimate, ...)
```

## Arguments

- x:

  A `loo` object, for example one returned by
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md),
  [`loo_subsample()`](https://mc-stan.org/loo/dev/reference/loo_subsample.md),
  [`loo_approximate_posterior()`](https://mc-stan.org/loo/dev/reference/loo_approximate_posterior.md),
  [`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md),
  etc.

- estimate:

  Which pointwise estimate to return. By default all are returned. The
  objects returned by the different functions
  ([`loo()`](https://mc-stan.org/loo/dev/reference/loo.md),
  [`loo_subsample()`](https://mc-stan.org/loo/dev/reference/loo_subsample.md),
  etc.) have slightly different estimates available. Typically at a
  minimum the estimates `elpd_loo`, `looic`, `mcse_elpd_loo`, `p_loo`,
  and `influence_pareto_k` will be available, but there may be others.

- ...:

  Currently ignored.

## Value

A vector of length equal to the number of observations.

## Examples

``` r
x <- loo(example_loglik_array())
pointwise(x, "elpd_loo")
#>  [1] -2.372054 -2.132756 -2.337100 -2.177556 -2.088849 -2.112027 -2.920085
#>  [8] -3.038168 -2.385516 -2.091878 -2.150381 -2.138015 -2.090231 -2.281848
#> [15] -2.241406 -2.448014 -4.517458 -4.989572 -2.330323 -4.711113 -2.466241
#> [22] -2.562420 -2.840067 -2.749268 -2.420624 -2.129352 -2.112648 -2.240625
#> [29] -3.220368 -2.505142 -2.642838 -2.145284
```
