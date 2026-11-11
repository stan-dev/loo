# Information Criteria (`ic`)

Computes the information criteria as -2 x log predictive density
(lppd_i) values. Inputs follow the same conventions as
[`measure_elpd()`](https://mc-stan.org/loo/pr/363/reference/measure_elpd.md)
and
[`measure_mlpd()`](https://mc-stan.org/loo/pr/363/reference/measure_mlpd.md).

## Usage

``` r
measure_ic(ylp, log_weights = NULL, pointwise = NULL)
```

## Arguments

- ylp:

  A numeric matrix or three-dimensional array of log predictive
  densities or probabilities:

  - **Matrix** (\\S \times n\\): \\S\\ posterior draws (chains merged)
    by \\n\\ observations.

  - **Array** (\\I \times C \times n\\): \\I\\ MCMC iterations per
    chain, \\C\\ chains, and \\n\\ observations. Converted to an \\S
    \times n\\ matrix internally.

- log_weights:

  Optional numeric matrix of unnormalized log-importance weights with
  dimensions \\S \times n\\. Weights are column-normalized before
  computing each per-observation contribution.

- pointwise:

  Optional numeric vector of precomputed pointwise contributions
  \\\mathrm{ic}\_i = -2 \cdot \mathrm{lppd}\_i\\. If provided, `ylp` and
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

## Examples

``` r
ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
measure_ic(ylp)
#> 
#> Computed from 2 draws by 2 observations.
#> 
#>    Estimate  SE
#> ic      3.6 1.2
```
