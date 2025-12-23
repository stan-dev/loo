# The number of posterior draws in a draws object.

The number of posterior draws in a draws object.

## Usage

``` r
.ndraws(x)

# S3 method for class 'matrix'
.ndraws(x)

# Default S3 method
.ndraws(x)
```

## Arguments

- x:

  A draws object with posterior draws.

## Value

An integer with the number of draws.

## Details

This is a generic function to return the total number of draws from an
arbitrary draws objects. The function is internal and should only be
used by developers to enable
[`loo_subsample()`](https://mc-stan.org/loo/reference/loo_subsample.md)
for arbitrary draws objects.
