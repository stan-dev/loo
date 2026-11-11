# Expected log pointwise predictive density (`elpd`)

Computes the expected log pointwise predictive density (ELPD) as the sum
of pointwise log predictive density contributions
(\\\mathrm{lppd}\_i\\), using
[`ptw_log_pred_density()`](https://mc-stan.org/loo/pr/363/reference/ptw_log_pred_density.md).
ELPD is returned on the utility scale (higher is better), consistent
with the sign convention used throughout this package.

## Usage

``` r
measure_elpd(ylp, log_weights = NULL, pointwise = NULL)
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

  Optional numeric vector of precomputed \\\mathrm{lppd}\_i\\ values.
  When supplied, `ylp` and `log_weights` are ignored.

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

\$\$\mathrm{elpd} = \sum\_{i=1}^{n} \mathrm{lppd}\_i,\$\$ where each
\\\mathrm{lppd}\_i\\ is computed by
[`ptw_log_pred_density()`](https://mc-stan.org/loo/pr/363/reference/ptw_log_pred_density.md).
The standard error is \\\sqrt{n}\\\mathrm{sd}(\mathrm{lppd}\_i)\\.

## See also

[`ptw_log_pred_density()`](https://mc-stan.org/loo/pr/363/reference/ptw_log_pred_density.md),
[`measure_mlpd()`](https://mc-stan.org/loo/pr/363/reference/measure_mlpd.md),
[`measure_ic()`](https://mc-stan.org/loo/pr/363/reference/measure_ic.md)

## Examples

``` r
ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
measure_elpd(ylp)
#> 
#> Computed from 2 draws by 2 observations.
#> 
#>      Estimate   SE
#> elpd     -1.8  0.6

# With unnormalized importance weights (e.g., PSIS-LOO)
lw <- matrix(log(c(0.7, 0.3, 0.6, 0.4)), nrow = 2)
measure_elpd(ylp, log_weights = lw)
#> 
#> Computed from 2 draws by 2 observations.
#> 
#>      Estimate   SE
#> elpd     -2.0  0.7

# From a draws x chains x observations array
LLarr <- example_loglik_array()
measure_elpd(LLarr)
#> 
#> Computed from 1000 draws by 32 observations.
#> 
#>      Estimate    SE
#> elpd    -80.3   3.2
```
