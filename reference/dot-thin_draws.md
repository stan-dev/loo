# Thin a draws object

Thin a draws object

## Usage

``` r
.thin_draws(draws, loo_approximation_draws)

# S3 method for class 'matrix'
.thin_draws(draws, loo_approximation_draws)

# S3 method for class 'numeric'
.thin_draws(draws, loo_approximation_draws)

# Default S3 method
.thin_draws(draws, loo_approximation_draws)
```

## Arguments

- draws:

  A draws object with posterior draws.

- loo_approximation_draws:

  The number of posterior draws to return (ie after thinning).

## Value

A thinned draws object.

## Details

This is a generic function to thin draws from arbitrary draws objects.
The function is internal and should only be used by developers to enable
[`loo_subsample()`](https://mc-stan.org/loo/reference/loo_subsample.md)
for arbitrary draws objects.
