# Predictive R-squared (`r2`)

Computes predictive R-squared as one minus the ratio of prediction MSE
to the empirical variance of `y`. The standard error is computed with a
first-order delta-method approximation.

## Usage

``` r
measure_r2(y, mupred, log_weights = NULL, pointwise = NULL)
```

## Arguments

- y:

  A numeric vector of observed outcomes.

- mupred:

  A numeric matrix of posterior expected predictions with dimensions
  draws x observations. A length-`n` vector is also accepted.

- log_weights:

  Optional numeric matrix of unnormalized log-importance weights with
  dimensions \\S \times n\\. Weights are column-normalized before
  computing each per-observation contribution.

- pointwise:

  Optional numeric vector of precomputed pointwise squared errors. If
  provided, `mupred` and `log_weights` are ignored. `y` is always
  required, because `r2` needs the variance of `y` as its baseline.

## Value

An object of class `"measure"`: a list with:

- `estimates`:

  A 1 \\\times\\ 2 matrix. The row is named after the measure, the
  columns are `Estimate` and `SE` (standard error).

- `pointwise`:

  An \\n \times 1\\ matrix of per-observation values. The column is
  named after the measure.

- `extra`:

  Optional. Auxiliary data for the standard error of a difference. Most
  measures do not set it.

Attributes `measure` (i.e., measure name) and `dims` (draws \\\times\\
observations) are also set. Use
[`print()`](https://rdrr.io/r/base/print.html) for a readable summary.

## Examples

``` r
y <- c(1, 2, 3)
mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
measure_r2(y, mupred)
#> 
#> Computed from 2 draws by 3 observations.
#> 
#>    Estimate    SE
#> r2    0.750 0.000
```
