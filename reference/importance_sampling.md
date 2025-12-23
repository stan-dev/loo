# A parent class for different importance sampling methods.

A parent class for different importance sampling methods.

## Usage

``` r
importance_sampling(log_ratios, method, ...)

# S3 method for class 'array'
importance_sampling(
  log_ratios,
  method,
  ...,
  r_eff = 1,
  cores = getOption("mc.cores", 1)
)

# S3 method for class 'matrix'
importance_sampling(
  log_ratios,
  method,
  ...,
  r_eff = 1,
  cores = getOption("mc.cores", 1)
)

# Default S3 method
importance_sampling(log_ratios, method, ..., r_eff = 1)
```

## Arguments

- log_ratios:

  An array, matrix, or vector of importance ratios on the log scale (for
  PSIS-LOO these are *negative* log-likelihood values). See the
  **Methods (by class)** section below for a detailed description of how
  to specify the inputs for each method.

- method:

  The importance sampling method to use. The following methods are
  implemented:

  - [`"psis"`](https://mc-stan.org/loo/reference/psis.md):
    Pareto-Smoothed Importance Sampling (PSIS). Default method.

  - [`"tis"`](https://mc-stan.org/loo/reference/tis.md): Truncated
    Importance Sampling (TIS) with truncation at `sqrt(S)`, where `S` is
    the number of posterior draws.

  - [`"sis"`](https://mc-stan.org/loo/reference/sis.md): Standard
    Importance Sampling (SIS).

- ...:

  Arguments passed on to the various methods.

- r_eff:

  Vector of relative effective sample size estimates containing one
  element per observation. The values provided should be the relative
  effective sample sizes of `1/exp(log_ratios)` (i.e., `1/ratios`). This
  is related to the relative efficiency of estimating the normalizing
  term in self-normalizing importance sampling. If `r_eff` is not
  provided then the reported PSIS effective sample sizes and Monte Carlo
  error estimates can be over-optimistic. If the posterior draws are
  (near) independent then `r_eff=1` can be used. `r_eff` has to be a
  scalar (same value is used for all observations) or a vector with
  length equal to the number of observations. The default value is 1.
  See the
  [`relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.md)
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
