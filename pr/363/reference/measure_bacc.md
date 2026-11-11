# Balanced classification accuracy (`bacc`)

Computes balanced accuracy by averaging class-specific mean accuracy,
giving each observed class equal weight regardless of class frequency.

## Usage

``` r
measure_bacc(y, mupred, log_weights = NULL, pointwise = NULL)
```

## Arguments

- y:

  An integer vector of observed class labels. Code it 0 or 1 when
  `mupred` is a matrix, and 1 to K when `mupred` is a 3-D array.

- mupred:

  A numeric array of posterior predictive means. For binary outcomes use
  a draws x observations matrix. For multiclass outcomes use a draws x
  observations x categories array.

- log_weights:

  Optional numeric matrix of unnormalized log-importance weights with
  dimensions \\S \times n\\. Weights are column-normalized before
  computing each per-observation contribution.

- pointwise:

  Optional numeric vector of precomputed pointwise contributions \\acc_i
  / (K n\_{y_i})\\, as stored by `measure_bacc()`. If provided, `mupred`
  and `log_weights` are ignored. `y` is always required, because it
  defines the classes.

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
y <- c(1L, 1L, 2L, 2L)
mupred <- array(
  c(0.8, 0.2, 0.7, 0.3, 0.3, 0.7, 0.2, 0.8),
  dim = c(1, 4, 2)
)
measure_bacc(y, mupred)
#> 
#> Computed from 1 draws by 4 observations.
#> 
#>      Estimate    SE
#> bacc    0.500 0.250
```
