# Generic function for K-fold cross-validation for developers

For developers of Bayesian modeling packages, **loo** includes a generic
function `kfold()` so that methods may be defined for K-fold CV without
name conflicts between packages. See, for example, the `kfold()` methods
in the **rstanarm** and **brms** packages.

The **Value** section below describes the objects that `kfold()` methods
should return in order to be compatible with
[`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md)
and the **loo** package print methods.

## Usage

``` r
kfold(x, ...)

is.kfold(x)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Arguments to pass to specific methods.

## Value

For developers defining a `kfold()` method for a class `"foo"`, the
`kfold.foo()` function should return a list with class
`c("kfold", "loo")` with at least the following named elements:

- `"estimates"`: A `1x2` matrix containing the ELPD estimate and its
  standard error. The matrix must have row name "`elpd_kfold`" and
  column names `"Estimate"` and `"SE"`.

- `"pointwise"`: A `Nx1` matrix with column name `"elpd_kfold"`
  containing the pointwise contributions for each data point.

It is important for the object to have at least these classes and
components so that it is compatible with other functions like
[`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md)
and [`print()`](https://rdrr.io/r/base/print.html) methods.
