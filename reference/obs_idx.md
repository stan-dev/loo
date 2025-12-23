# Get observation indices used in subsampling

Get observation indices used in subsampling

## Usage

``` r
obs_idx(x, rep = TRUE)
```

## Arguments

- x:

  A `psis_loo_ss` object.

- rep:

  If sampling with replacement is used, an observation can have multiple
  samples and these are then repeated in the returned object if
  `rep=TRUE` (e.g., a vector `c(1,1,2)` indicates that observation 1 has
  been subampled two times). If `rep=FALSE` only the unique indices are
  returned.

## Value

An integer vector.
