# Extract importance sampling weights

Extract importance sampling weights

## Usage

``` r
# S3 method for class 'importance_sampling'
weights(object, ..., log = TRUE, normalize = TRUE)
```

## Arguments

- object:

  An object returned by
  [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md),
  [`tis()`](https://mc-stan.org/loo/dev/reference/tis.md), or
  [`sis()`](https://mc-stan.org/loo/dev/reference/sis.md).

- ...:

  Ignored.

- log:

  Should the weights be returned on the log scale? Defaults to `TRUE`.

- normalize:

  Should the weights be normalized? Defaults to `TRUE`.

## Value

The [`weights()`](https://rdrr.io/r/stats/weights.html) method returns
an object with the same dimensions as the `log_weights` component of
`object`. The `normalize` and `log` arguments control whether the
returned weights are normalized and whether or not to return them on the
log scale.

## Examples

``` r
# See the examples at help("psis")
```
