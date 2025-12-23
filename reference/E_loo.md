# Compute weighted expectations

The `E_loo()` function computes weighted expectations (means, variances,
quantiles) using the importance weights obtained from the
[PSIS](https://mc-stan.org/loo/reference/psis.md) smoothing procedure.
The expectations estimated by the `E_loo()` function assume that the
PSIS approximation is working well. **A small [Pareto
k](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md) estimate is
necessary, but not sufficient, for `E_loo()` to give reliable
estimates**. If the `log_ratios` argument is provided, `E_loo()` also
computes a function specific Pareto k diagnostic, which must also be
small for a reliable estimate. See more details below.

## Usage

``` r
E_loo(x, psis_object, ...)

# Default S3 method
E_loo(
  x,
  psis_object,
  ...,
  type = c("mean", "variance", "sd", "quantile"),
  probs = NULL,
  log_ratios = NULL
)

# S3 method for class 'matrix'
E_loo(
  x,
  psis_object,
  ...,
  type = c("mean", "variance", "sd", "quantile"),
  probs = NULL,
  log_ratios = NULL
)
```

## Arguments

- x:

  A numeric vector or matrix.

- psis_object:

  An object returned by
  [`psis()`](https://mc-stan.org/loo/reference/psis.md).

- ...:

  Arguments passed to individual methods.

- type:

  The type of expectation to compute. The options are `"mean"`,
  `"variance"`, `"sd"`, and `"quantile"`.

- probs:

  For computing quantiles, a vector of probabilities.

- log_ratios:

  Optionally, a vector or matrix (the same dimensions as `x`) of raw
  (not smoothed) log ratios. If working with log-likelihood values, the
  log ratios are the **negative** of those values. If `log_ratios` is
  specified we are able to compute more accurate [Pareto
  k](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  diagnostics specific to `E_loo()`.

## Value

A named list with the following components:

- `value`:

  The result of the computation.

  For the matrix method, `value` is a vector with `ncol(x)` elements,
  with one exception: when `type="quantile"` and multiple values are
  specified in `probs` the `value` component of the returned object is a
  `length(probs)` by `ncol(x)` matrix.

  For the default/vector method the `value` component is scalar, with
  one exception: when `type="quantile"` and multiple values are
  specified in `probs` the `value` component is a vector with
  `length(probs)` elements.

- `pareto_k`:

  Function-specific diagnostic.

  For the matrix method it will be a vector of length `ncol(x)`
  containing estimates of the shape parameter \\k\\ of the generalized
  Pareto distribution. For the default/vector method, the estimate is a
  scalar. If `log_ratios` is not specified when calling `E_loo()`, the
  smoothed log-weights are used to estimate Pareto-k's, which may
  produce optimistic estimates.

  For `type="mean"`, `type="var"`, and `type="sd"`, the returned
  Pareto-k is usually the maximum of the Pareto-k's for the left and
  right tail of \\hr\\ and the right tail of \\r\\, where \\r\\ is the
  importance ratio and \\h=x\\ for `type="mean"` and \\h=x^2\\ for
  `type="var"` and `type="sd"`. If \\h\\ is binary, constant, or not
  finite, or if `type="quantile"`, the returned Pareto-k is the Pareto-k
  for the right tail of \\r\\.

## Examples

``` r
# \donttest{
if (requireNamespace("rstanarm", quietly = TRUE)) {
# Use rstanarm package to quickly fit a model and get both a log-likelihood
# matrix and draws from the posterior predictive distribution
library("rstanarm")

# data from help("lm")
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
d <- data.frame(
  weight = c(ctl, trt),
  group = gl(2, 10, 20, labels = c("Ctl","Trt"))
)
fit <- stan_glm(weight ~ group, data = d, refresh = 0)
yrep <- posterior_predict(fit)
dim(yrep)

log_ratios <- -1 * log_lik(fit)
dim(log_ratios)

r_eff <- relative_eff(exp(-log_ratios), chain_id = rep(1:4, each = 1000))
psis_object <- psis(log_ratios, r_eff = r_eff, cores = 2)

E_loo(yrep, psis_object, type = "mean")
E_loo(yrep, psis_object, type = "var")
E_loo(yrep, psis_object, type = "sd")
E_loo(yrep, psis_object, type = "quantile", probs = 0.5) # median
E_loo(yrep, psis_object, type = "quantile", probs = c(0.1, 0.9))

# We can get more accurate Pareto k diagnostic if we also provide
# the log_ratios argument
E_loo(yrep, psis_object, type = "mean", log_ratios = log_ratios)
}
#> Loading required package: Rcpp
#> This is rstanarm version 2.32.2
#> - See https://mc-stan.org/rstanarm/articles/priors for changes to default priors!
#> - Default priors may change, so it's safest to specify priors, even if equivalent to the defaults.
#> - For execution on a local, multicore CPU with excess RAM we recommend calling
#>   options(mc.cores = parallel::detectCores())
#> $value
#>  [1] 5.141494 4.952522 5.023677 4.911770 5.085666 5.064760 5.010091 5.098857
#>  [9] 5.004081 5.010727 4.623576 4.730398 4.689493 4.791352 4.535483 4.727869
#> [17] 4.533920 4.645670 4.700269 4.643684
#> 
#> $pareto_k
#>  [1] 0.29359847 0.15973717 0.12864056 0.32680903 0.23337975 0.27264338
#>  [7] 0.19035143 0.31173918 0.16756065 0.15895287 0.09213184 0.10872512
#> [13] 0.12425409 0.34636785 0.33291067 0.30275686 0.47487236 0.06150852
#> [19] 0.05578757 0.03591544
#> 
# }
```
