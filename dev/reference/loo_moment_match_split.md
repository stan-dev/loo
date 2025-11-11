# Split moment matching for efficient approximate leave-one-out cross-validation (LOO)

A function that computes the split moment matching importance sampling
loo. Takes in the moment matching total transformation, transforms only
half of the draws, and computes a single elpd using multiple importance
sampling.

## Usage

``` r
loo_moment_match_split(
  x,
  upars,
  cov,
  total_shift,
  total_scaling,
  total_mapping,
  i,
  log_prob_upars,
  log_lik_i_upars,
  r_eff_i,
  cores,
  is_method,
  ...
)
```

## Arguments

- x:

  A fitted model object.

- upars:

  A matrix containing the model parameters in unconstrained space where
  they can have any real value.

- cov:

  Logical; Indicate whether to match the covariance matrix of the
  samples or not. If `FALSE`, only the mean and marginal variances are
  matched.

- total_shift:

  A vector representing the total shift made by the moment matching
  algorithm.

- total_scaling:

  A vector representing the total scaling of marginal variance made by
  the moment matching algorithm.

- total_mapping:

  A vector representing the total covariance transformation made by the
  moment matching algorithm.

- i:

  Observation index.

- log_prob_upars:

  A function that takes arguments `x` and `upars` and returns a matrix
  of log-posterior density values of the unconstrained posterior draws
  passed via `upars`.

- log_lik_i_upars:

  A function that takes arguments `x`, `upars`, and `i` and returns a
  vector of log-likeliood draws of the `i`th observation based on the
  unconstrained posterior draws passed via `upars`.

- r_eff_i:

  MCMC relative effective sample size of the `i`'th log likelihood
  draws.

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

- is_method:

  The importance sampling method to use. The following methods are
  implemented:

  - [`"psis"`](https://mc-stan.org/loo/dev/reference/psis.md):
    Pareto-Smoothed Importance Sampling (PSIS). Default method.

  - [`"tis"`](https://mc-stan.org/loo/dev/reference/tis.md): Truncated
    Importance Sampling (TIS) with truncation at `sqrt(S)`, where `S` is
    the number of posterior draws.

  - [`"sis"`](https://mc-stan.org/loo/dev/reference/sis.md): Standard
    Importance Sampling (SIS).

- ...:

  Further arguments passed to the custom functions documented above.

## Value

A list containing the updated log-importance weights and log-likelihood
values. Also returns the updated MCMC effective sample size and the
integrand-specific log-importance weights.

## References

Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2021).
Implicitly adaptive importance sampling. *Statistics and Computing*, 31,
16. doi:10.1007/s11222-020-09982-2. arXiv preprint arXiv:1906.08850.

## See also

[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md),
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
