# Pointwise log predictive density (`lppd_i`)

Computes pointwise log predictive density contributions from a matrix of
log predictive densities/probabilities for posterior draws. When PSIS
log-weights are supplied, they are used to form a weighted log-sum-exp
per observation.

## Usage

``` r
ptw_log_pred_density(ylp, psis_log_weights = NULL)
```

## Arguments

- ylp:

  A numeric matrix of log predictive densities/probabilities with
  dimensions draws x observations.

- psis_log_weights:

  Optional numeric matrix of normalized PSIS log weights with the same
  dimensions as `ylp`. Each column must sum to 1 on the probability
  scale.

## Value

A numeric vector of length `ncol(ylp)` with pointwise log predictive
density values.

## Examples

``` r
ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
ptw_log_pred_density(ylp)
#> [1] -1.203973 -0.597837

lw <- matrix(log(c(0.7, 0.3, 0.6, 0.4)), nrow = 2)
ptw_log_pred_density(ylp, lw)
#> [1] -1.3470736 -0.6931472
```
