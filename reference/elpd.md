# Generic (expected) log-predictive density

The `elpd()` methods for arrays and matrices can compute the expected
log pointwise predictive density for a new dataset or the log pointwise
predictive density of the observed data (an overestimate of the elpd).

## Usage

``` r
elpd(x, ...)

# S3 method for class 'array'
elpd(x, ...)

# S3 method for class 'matrix'
elpd(x, ...)
```

## Arguments

- x:

  A log-likelihood array or matrix. The **Methods (by class)** section,
  below, has detailed descriptions of how to specify the inputs for each
  method.

- ...:

  Currently ignored.

## Details

The `elpd()` function is an S3 generic and methods are provided for 3-D
pointwise log-likelihood arrays and matrices.

## Methods (by class)

- `elpd(array)`: An \\I\\ by \\C\\ by \\N\\ array, where \\I\\ is the
  number of MCMC iterations per chain, \\C\\ is the number of chains,
  and \\N\\ is the number of data points.

- `elpd(matrix)`: An \\S\\ by \\N\\ matrix, where \\S\\ is the size of
  the posterior sample (with all chains merged) and \\N\\ is the number
  of data points.

## See also

The vignette *Holdout validation and K-fold cross-validation of Stan
programs with the loo package* for demonstrations of using the `elpd()`
methods.

## Examples

``` r
# Calculate the lpd of the observed data
LLarr <- example_loglik_array()
elpd(LLarr)
#> 
#> Computed from 1000 by 32 log-likelihood matrix using the generic elpd function
#> 
#>      Estimate  SE
#> elpd    -80.3 3.2
#> ic      160.5 6.5
```
