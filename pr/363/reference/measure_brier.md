# Brier score (`brier`)

Computes the Brier score for binary outcomes as squared error between
the observed label and predicted event probability.

## Usage

``` r
measure_brier(y, ypred, log_weights = NULL, pointwise = NULL)
```

## Arguments

- y:

  A numeric vector of binary outcomes coded as 0 or 1.

- ypred:

  A numeric matrix of posterior predictive probabilities with dimensions
  draws x observations.

- log_weights:

  Optional numeric matrix of unnormalized log-importance weights with
  dimensions \\S \times n\\. Weights are column-normalized before
  computing each per-observation contribution.

- pointwise:

  Optional numeric vector of precomputed pointwise Brier scores. If
  provided, `y`, `ypred`, and `log_weights` are ignored.

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
y <- c(1, 0, 1)
ypred <- matrix(c(0.8, 0.2, 0.7, 0.9, 0.4, 0.6), nrow = 2)
measure_brier(y, ypred)
#> 
#> Computed from 2 draws by 3 observations.
#> 
#>       Estimate    SE
#> brier    0.380 0.130
```
