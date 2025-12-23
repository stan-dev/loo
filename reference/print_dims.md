# Print dimensions of log-likelihood or log-weights matrix

Print dimensions of log-likelihood or log-weights matrix

## Usage

``` r
print_dims(x, ...)

# S3 method for class 'importance_sampling'
print_dims(x, ...)

# S3 method for class 'psis_loo'
print_dims(x, ...)

# S3 method for class 'importance_sampling_loo'
print_dims(x, ...)

# S3 method for class 'waic'
print_dims(x, ...)

# S3 method for class 'kfold'
print_dims(x, ...)

# S3 method for class 'psis_loo_ss'
print_dims(x, ...)
```

## Arguments

- x:

  The object returned by
  [`psis()`](https://mc-stan.org/loo/reference/psis.md),
  [`loo()`](https://mc-stan.org/loo/reference/loo.md), or
  [`waic()`](https://mc-stan.org/loo/reference/waic.md).

- ...:

  Ignored.
