# Scaled Ranked Probability Score (SRPS, SCRPS)

A convenience wrapper around
[`measure_rps()`](https://mc-stan.org/loo/pr/363/reference/measure_rps.md)
with `scaled = TRUE`. Computes the scaled ranked probability score
(SRPS) for discrete outcomes or the scaled continuously ranked
probability score (SCRPS) for continuous outcomes. Scaling makes the
score invariant to the spread of the predictive distribution, which can
be useful when comparing models across observations with very different
predictive uncertainties.

## Usage

``` r
measure_srps(y, ypred, log_weights = NULL, pointwise = NULL)
```

## Arguments

- y:

  A numeric vector of \\n\\ observed outcomes. May be integer-valued
  (for RPS/SRPS) or continuous (for CRPS/SCRPS).

- ypred:

  A numeric matrix of posterior predictive draws with dimensions \\S
  \times n\\ (draws × observations).

- log_weights:

  Optional numeric matrix of unnormalized log-importance weights with
  dimensions \\S \times n\\. Weights are column-normalized before
  computing each per-observation contribution.

- pointwise:

  Optional numeric vector of precomputed pointwise values, on the scale
  of the measure that is returned: RPS/CRPS when `scaled = FALSE` and
  SRPS/SCRPS when `scaled = TRUE`. The values are used as they are, with
  no further transformation; in particular, unscaled values are not
  converted to the scaled variant. If provided, `y`, `ypred`, and
  `log_weights` are ignored.

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

## Details

See
[`measure_rps()`](https://mc-stan.org/loo/pr/363/reference/measure_rps.md)
for full details on arguments, estimation, and references.

## Examples

``` r
y <- c(2L, 1L, 3L)
ypred <- matrix(c(2, 1, 2, 3, 1, 3), nrow = 2)
measure_srps(y, ypred)
#> 
#> Computed from 2 draws by 3 observations.
#> 
#>      Estimate    SE
#> srps    -0.95  0.29
```
