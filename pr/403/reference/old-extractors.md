# Extractor methods

These are only defined in order to deprecate with a warning (rather than
remove and break backwards compatibility) the old way of accessing the
point estimates in a `"psis_loo"` or `"psis"` object. The new way as of
v2.0.0 is to get them from the `"estimates"` component of the object.

## Usage

``` r
# S3 method for class 'loo'
x[i]

# S3 method for class 'loo'
x[[i, exact = TRUE]]

# S3 method for class 'loo'
x$name
```

## Arguments

- x, i, exact, name:

  See [Extract](https://rdrr.io/r/base/Extract.html).
