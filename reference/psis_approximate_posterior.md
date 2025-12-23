# Diagnostics for Laplace and ADVI approximations and Laplace-loo and ADVI-loo

Diagnostics for Laplace and ADVI approximations and Laplace-loo and
ADVI-loo

## Usage

``` r
psis_approximate_posterior(
  log_p = NULL,
  log_g = NULL,
  log_liks = NULL,
  cores,
  save_psis,
  ...,
  log_q = NULL
)
```

## Arguments

- log_p:

  The log-posterior (target) evaluated at S samples from the proposal
  distribution (g). A vector of length S.

- log_g:

  The log-density (proposal) evaluated at S samples from the proposal
  distribution (g). A vector of length S.

- log_liks:

  A log-likelihood matrix of size S \* N, where N is the number of
  observations and S is the number of samples from q. See
  [`loo.matrix()`](https://mc-stan.org/loo/reference/loo.md) for
  details. Default is `NULL`. Then only the posterior is evaluated using
  the k_hat diagnostic.

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

- save_psis:

  Should the `psis` object created internally by
  [`loo()`](https://mc-stan.org/loo/reference/loo.md) be saved in the
  returned object? The
  [`loo()`](https://mc-stan.org/loo/reference/loo.md) function calls
  [`psis()`](https://mc-stan.org/loo/reference/psis.md) internally but
  by default discards the (potentially large) `psis` object after using
  it to compute the LOO-CV summaries. Setting `save_psis=TRUE` will add
  a `psis_object` component to the list returned by `loo`. This is
  useful if you plan to use the
  [`E_loo()`](https://mc-stan.org/loo/reference/E_loo.md) function to
  compute weighted expectations after running `loo`. Several functions
  in the bayesplot package also accept `psis` objects.

- log_q:

  Deprecated argument name (the same as log_g).

## Value

If log likelihoods are supplied, the function returns a `"loo"` object,
otherwise the function returns a `"psis"` object.

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

[`loo()`](https://mc-stan.org/loo/reference/loo.md) and
[`psis()`](https://mc-stan.org/loo/reference/psis.md)
