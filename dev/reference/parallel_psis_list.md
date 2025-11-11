# Parallel psis list computations

Parallel psis list computations

## Usage

``` r
parallel_psis_list(
  N,
  .loo_i,
  .llfun,
  data,
  draws,
  r_eff,
  save_psis,
  cores,
  ...
)

parallel_importance_sampling_list(
  N,
  .loo_i,
  .llfun,
  data,
  draws,
  r_eff,
  save_psis,
  cores,
  method,
  ...
)
```

## Arguments

- N:

  The total number of observations (i.e. `nrow(data)`).

- .loo_i:

  The function used to compute individual loo contributions.

- .llfun:

  See `llfun` in
  [`loo.function()`](https://mc-stan.org/loo/dev/reference/loo.md).

- data, draws, ...:

  For the
  [`loo.function()`](https://mc-stan.org/loo/dev/reference/loo.md)
  method and the
  [`loo_i()`](https://mc-stan.org/loo/dev/reference/loo.md) function,
  these are the data, posterior draws, and other arguments to pass to
  the log-likelihood function. See the **Methods (by class)** section
  below for details on how to specify these arguments.

- r_eff:

  Vector of relative effective sample size estimates for the likelihood
  (`exp(log_lik)`) of each observation. This is related to the relative
  efficiency of estimating the normalizing term in self-normalized
  importance sampling when using posterior draws obtained with MCMC. If
  MCMC draws are used and `r_eff` is not provided then the reported PSIS
  effective sample sizes and Monte Carlo error estimates can be
  over-optimistic. If the posterior draws are (near) independent then
  `r_eff=1` can be used. `r_eff` has to be a scalar (same value is used
  for all observations) or a vector with length equal to the number of
  observations. The default value is 1. See the
  [`relative_eff()`](https://mc-stan.org/loo/dev/reference/relative_eff.md)
  helper functions for help computing `r_eff`.

- save_psis:

  Should the `psis` object created internally by
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) be saved in
  the returned object? The
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) function calls
  [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md) internally
  but by default discards the (potentially large) `psis` object after
  using it to compute the LOO-CV summaries. Setting `save_psis=TRUE`
  will add a `psis_object` component to the list returned by `loo`. This
  is useful if you plan to use the
  [`E_loo()`](https://mc-stan.org/loo/dev/reference/E_loo.md) function
  to compute weighted expectations after running `loo`. Several
  functions in the bayesplot package also accept `psis` objects.

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

- method:

  See `is_method` for
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md)

## Details

Refactored function to handle parallel computations for psis_list
