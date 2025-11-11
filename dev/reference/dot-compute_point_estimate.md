# Compute a point estimate from a draws object

Compute a point estimate from a draws object

## Usage

``` r
.compute_point_estimate(draws)

# S3 method for class 'matrix'
.compute_point_estimate(draws)

# Default S3 method
.compute_point_estimate(draws)
```

## Arguments

- draws:

  A draws object with draws from the posterior.

## Value

A 1 by P matrix with point estimates from a draws object.

## Details

This is a generic function to compute point estimates from draws
objects. The function is internal and should only be used by developers
to enable
[`loo_subsample()`](https://mc-stan.org/loo/dev/reference/loo_subsample.md)
for arbitrary draws objects.
