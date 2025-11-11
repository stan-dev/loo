# Truncated importance sampling (TIS)

Implementation of truncated (self-normalized) importance sampling (TIS),
truncated at S^(1/2) as recommended by Ionides (2008).

## Usage

``` r
tis(log_ratios, ...)

# S3 method for class 'array'
tis(log_ratios, ..., r_eff = 1, cores = getOption("mc.cores", 1))

# S3 method for class 'matrix'
tis(log_ratios, ..., r_eff = 1, cores = getOption("mc.cores", 1))

# Default S3 method
tis(log_ratios, ..., r_eff = 1)
```

## Arguments

- log_ratios:

  An array, matrix, or vector of importance ratios on the log scale (for
  Importance sampling LOO, these are *negative* log-likelihood values).
  See the **Methods (by class)** section below for a detailed
  description of how to specify the inputs for each method.

- ...:

  Arguments passed on to the various methods.

- r_eff:

  Vector of relative effective sample size estimates containing one
  element per observation. The values provided should be the relative
  effective sample sizes of `1/exp(log_ratios)` (i.e., `1/ratios`). This
  is related to the relative efficiency of estimating the normalizing
  term in self-normalizing importance sampling. If `r_eff` is not
  provided then the reported (T)IS effective sample sizes and Monte
  Carlo error estimates can be over-optimistic. If the posterior draws
  are (near) independent then `r_eff=1` can be used. `r_eff` has to be a
  scalar (same value is used for all observations) or a vector with
  length equal to the number of observations. The default value is 1.
  See the
  [`relative_eff()`](https://mc-stan.org/loo/dev/reference/relative_eff.md)
  helper function for computing `r_eff`.

- cores:

  The number of cores to use for parallelization. This defaults to the
  option `mc.cores` which can be set for an entire R session by
  `options(mc.cores = NUMBER)`. The old option `loo.cores` is now
  deprecated but will be given precedence over `mc.cores` until
  `loo.cores` is removed in a future release. **As of version 2.0.0 the
  default is now 1 core if `mc.cores` is not set**, but we recommend
  using as many (or close to as many) cores as possible.

  - Note for Windows 10 users: it is **strongly**
    [recommended](https://github.com/stan-dev/loo/issues/94) to avoid
    using the `.Rprofile` file to set `mc.cores` (using the `cores`
    argument or setting `mc.cores` interactively or in a script is
    fine).

## Value

The `tis()` methods return an object of class `"tis"`, which is a named
list with the following components:

- `log_weights`:

  Vector or matrix of smoothed (and truncated) but *unnormalized* log
  weights. To get normalized weights use the
  [`weights()`](https://mc-stan.org/loo/dev/reference/weights.importance_sampling.md)
  method provided for objects of class `tis`.

- `diagnostics`:

  A named list containing one vector:

  - `pareto_k`: Not used in `tis`, all set to 0.

  - `n_eff`: Effective sample size estimates.

Objects of class `"tis"` also have the following
[attributes](https://rdrr.io/r/base/attributes.html):

- `norm_const_log`:

  Vector of precomputed values of `colLogSumExps(log_weights)` that are
  used internally by the
  [`weights()`](https://rdrr.io/r/stats/weights.html)method to normalize
  the log weights.

- `r_eff`:

  If specified, the user's `r_eff` argument.

- `tail_len`:

  Not used for `tis`.

- `dims`:

  Integer vector of length 2 containing `S` (posterior sample size) and
  `N` (number of observations).

- `method`:

  Method used for importance sampling, here `tis`.

## Methods (by class)

- `tis(array)`: An \\I\\ by \\C\\ by \\N\\ array, where \\I\\ is the
  number of MCMC iterations per chain, \\C\\ is the number of chains,
  and \\N\\ is the number of data points.

- `tis(matrix)`: An \\S\\ by \\N\\ matrix, where \\S\\ is the size of
  the posterior sample (with all chains merged) and \\N\\ is the number
  of data points.

- `tis(default)`: A vector of length \\S\\ (posterior sample size).

## References

Ionides, Edward L. (2008). Truncated importance sampling. *Journal of
Computational and Graphical Statistics* 17(2): 295–311.

## See also

- [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md) for
  approximate LOO-CV using PSIS.

- [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) for
  approximate LOO-CV.

- [pareto-k-diagnostic](https://mc-stan.org/loo/dev/reference/pareto-k-diagnostic.md)
  for PSIS diagnostics.

## Examples

``` r
log_ratios <- -1 * example_loglik_array()
r_eff <- relative_eff(exp(-log_ratios))
tis_result <- tis(log_ratios, r_eff = r_eff)
str(tis_result)
#> List of 2
#>  $ log_weights: num [1:1000, 1:32] 2.37 2.12 2.24 2.41 2.25 ...
#>  $ diagnostics:List of 3
#>   ..$ pareto_k: num [1:32] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ n_eff   : num [1:32] 910 937 938 901 907 ...
#>   ..$ r_eff   : num [1:32] 0.942 0.954 0.977 0.919 0.923 ...
#>  - attr(*, "norm_const_log")= num [1:32] 9.28 9.04 9.24 9.09 9 ...
#>  - attr(*, "tail_len")= num [1:32] 98 98 96 99 99 101 99 100 102 98 ...
#>  - attr(*, "r_eff")= num [1:32] 0.942 0.954 0.977 0.919 0.923 ...
#>  - attr(*, "dims")= int [1:2] 1000 32
#>  - attr(*, "method")= chr "tis"
#>  - attr(*, "class")= chr [1:3] "tis" "importance_sampling" "list"

# extract smoothed weights
lw <- weights(tis_result) # default args are log=TRUE, normalize=TRUE
ulw <- weights(tis_result, normalize=FALSE) # unnormalized log-weights

w <- weights(tis_result, log=FALSE) # normalized weights (not log-weights)
uw <- weights(tis_result, log=FALSE, normalize = FALSE) # unnormalized weights
```
