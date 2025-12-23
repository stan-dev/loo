# Standard importance sampling (SIS)

Implementation of standard importance sampling (SIS).

## Usage

``` r
sis(log_ratios, ...)

# S3 method for class 'array'
sis(log_ratios, ..., r_eff = NULL, cores = getOption("mc.cores", 1))

# S3 method for class 'matrix'
sis(log_ratios, ..., r_eff = NULL, cores = getOption("mc.cores", 1))

# Default S3 method
sis(log_ratios, ..., r_eff = NULL)
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
  term in self-normalizing importance sampling. See the
  [`relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.md)
  helper function for computing `r_eff`. If using `psis` with draws of
  the `log_ratios` not obtained from MCMC then the warning message
  thrown when not specifying `r_eff` can be disabled by setting `r_eff`
  to `NA`.

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

The `sis()` methods return an object of class `"sis"`, which is a named
list with the following components:

- `log_weights`:

  Vector or matrix of smoothed but *unnormalized* log weights. To get
  normalized weights use the
  [`weights()`](https://mc-stan.org/loo/reference/weights.importance_sampling.md)
  method provided for objects of class `sis`.

- `diagnostics`:

  A named list containing one vector:

  - `pareto_k`: Not used in `sis`, all set to 0.

  - `n_eff`: effective sample size estimates.

Objects of class `"sis"` also have the following
[attributes](https://rdrr.io/r/base/attributes.html):

- `norm_const_log`:

  Vector of precomputed values of `colLogSumExps(log_weights)` that are
  used internally by the `weights` method to normalize the log weights.

- `r_eff`:

  If specified, the user's `r_eff` argument.

- `tail_len`:

  Not used for `sis`.

- `dims`:

  Integer vector of length 2 containing `S` (posterior sample size) and
  `N` (number of observations).

- `method`:

  Method used for importance sampling, here `sis`.

## Methods (by class)

- `sis(array)`: An \\I\\ by \\C\\ by \\N\\ array, where \\I\\ is the
  number of MCMC iterations per chain, \\C\\ is the number of chains,
  and \\N\\ is the number of data points.

- `sis(matrix)`: An \\S\\ by \\N\\ matrix, where \\S\\ is the size of
  the posterior sample (with all chains merged) and \\N\\ is the number
  of data points.

- `sis(default)`: A vector of length \\S\\ (posterior sample size).

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4
([journal
version](https://link.springer.com/article/10.1007/s11222-016-9696-4),
[preprint arXiv:1507.04544](https://arxiv.org/abs/1507.04544)).

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72):1-58. [PDF](https://jmlr.org/papers/v25/19-556.html)

## See also

- [`psis()`](https://mc-stan.org/loo/reference/psis.md) for approximate
  LOO-CV using PSIS.

- [`loo()`](https://mc-stan.org/loo/reference/loo.md) for approximate
  LOO-CV.

- [pareto-k-diagnostic](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  for PSIS diagnostics.

## Examples

``` r
log_ratios <- -1 * example_loglik_array()
r_eff <- relative_eff(exp(-log_ratios))
sis_result <- sis(log_ratios, r_eff = r_eff)
str(sis_result)
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
#>  - attr(*, "method")= chr "sis"
#>  - attr(*, "class")= chr [1:3] "sis" "importance_sampling" "list"

# extract smoothed weights
lw <- weights(sis_result) # default args are log=TRUE, normalize=TRUE
ulw <- weights(sis_result, normalize=FALSE) # unnormalized log-weights

w <- weights(sis_result, log=FALSE) # normalized weights (not log-weights)
uw <- weights(sis_result, log=FALSE, normalize = FALSE) # unnormalized weights
```
