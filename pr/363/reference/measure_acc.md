# Classification accuracy (`acc`)

Computes pointwise and average classification accuracy for binary or
multiclass outcomes. The class probabilities are first averaged over the
draws, using `log_weights` if given. For binary outcomes, the averaged
probability is thresholded at 0.5. For multiclass outcomes, the
predicted class is the most likely category of the averaged
probabilities, via
[`which.max()`](https://rdrr.io/r/base/which.min.html).

## Usage

``` r
measure_acc(y, mupred, log_weights = NULL, pointwise = NULL)
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

  Optional numeric vector of precomputed pointwise accuracy
  contributions. If provided, `y`, `mupred`, and `log_weights` are
  ignored.

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
y <- c(1L, 0L, 1L)
mupred <- matrix(c(0.8, 0.3, 0.7, 0.6, 0.4, 0.9), nrow = 2)
measure_acc(y, mupred)
#> 
#> Computed from 2 draws by 3 observations.
#> 
#>     Estimate    SE
#> acc    0.667 0.272
```
