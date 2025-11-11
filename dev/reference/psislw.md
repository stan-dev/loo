# Pareto smoothed importance sampling (deprecated, old version)

As of version `2.0.0` this function is **deprecated**. Please use the
[`psis()`](https://mc-stan.org/loo/dev/reference/psis.md) function for
the new PSIS algorithm.

## Usage

``` r
psislw(
  lw,
  wcp = 0.2,
  wtrunc = 3/4,
  cores = getOption("mc.cores", 1),
  llfun = NULL,
  llargs = NULL,
  ...
)
```

## Arguments

- lw:

  A matrix or vector of log weights. For computing LOO, `lw = -log_lik`,
  the *negative* of an \\S\\ (simulations) by \\N\\ (data points)
  pointwise log-likelihood matrix.

- wcp:

  The proportion of importance weights to use for the generalized Pareto
  fit. The `100*wcp`\\ from which to estimate the parameters of the
  generalized Pareto distribution.

- wtrunc:

  For truncating very large weights to \\S\\^`wtrunc`. Set to zero for
  no truncation.

- cores:

  The number of cores to use for parallelization. This defaults to the
  option `mc.cores` which can be set for an entire R session by
  `options(mc.cores = NUMBER)`, the old option `loo.cores` is now
  deprecated but will be given precedence over `mc.cores` until it is
  removed. **As of version 2.0.0, the default is now 1 core if
  `mc.cores` is not set, but we recommend using as many (or close to as
  many) cores as possible.**

- llfun, llargs:

  See [`loo.function()`](https://mc-stan.org/loo/dev/reference/loo.md).

- ...:

  Ignored when `psislw()` is called directly. The `...` is only used
  internally when `psislw()` is called by the
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) function.

## Value

A named list with components `lw_smooth` (modified log weights) and
`pareto_k` (estimated generalized Pareto shape parameter(s) k).

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

[pareto-k-diagnostic](https://mc-stan.org/loo/dev/reference/pareto-k-diagnostic.md)
for PSIS diagnostics.
