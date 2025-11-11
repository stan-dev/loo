# Moment matching for efficient approximate leave-one-out cross-validation (LOO)

Moment matching algorithm for updating a loo object when Pareto k
estimates are large.

## Usage

``` r
loo_moment_match(x, ...)

# Default S3 method
loo_moment_match(
  x,
  loo,
  post_draws,
  log_lik_i,
  unconstrain_pars,
  log_prob_upars,
  log_lik_i_upars,
  max_iters = 30L,
  k_threshold = NULL,
  split = TRUE,
  cov = TRUE,
  cores = getOption("mc.cores", 1),
  ...
)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Further arguments passed to the custom functions documented above.

- loo:

  A loo object to be modified.

- post_draws:

  A function the takes `x` as the first argument and returns a matrix of
  posterior draws of the model parameters.

- log_lik_i:

  A function that takes `x` and `i` and returns a matrix (one column per
  chain) or a vector (all chains stacked) of log-likelihood draws of the
  `i`th observation based on the model `x`. If the draws are obtained
  using MCMC, the matrix with MCMC chains separated is preferred.

- unconstrain_pars:

  A function that takes arguments `x`, and `pars` and returns posterior
  draws on the unconstrained space based on the posterior draws on the
  constrained space passed via `pars`.

- log_prob_upars:

  A function that takes arguments `x` and `upars` and returns a matrix
  of log-posterior density values of the unconstrained posterior draws
  passed via `upars`.

- log_lik_i_upars:

  A function that takes arguments `x`, `upars`, and `i` and returns a
  vector of log-likelihood draws of the `i`th observation based on the
  unconstrained posterior draws passed via `upars`.

- max_iters:

  Maximum number of moment matching iterations. Usually this does not
  need to be modified. If the maximum number of iterations is reached,
  there will be a warning, and increasing `max_iters` may improve
  accuracy.

- k_threshold:

  Threshold value for Pareto k values above which the moment matching
  algorithm is used. The default value is `min(1 - 1/log10(S), 0.7)`,
  where `S` is the sample size.

- split:

  Logical; Indicate whether to do the split transformation or not at the
  end of moment matching for each LOO fold.

- cov:

  Logical; Indicate whether to match the covariance matrix of the
  samples or not. If `FALSE`, only the mean and marginal variances are
  matched.

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

The `loo_moment_match()` methods return an updated `loo` object. The
structure of the updated `loo` object is similar, but the method also
stores the original Pareto k diagnostic values in the diagnostics field.

## Details

The `loo_moment_match()` function is an S3 generic and we provide a
default method that takes as arguments user-specified functions
`post_draws`, `log_lik_i`, `unconstrain_pars`, `log_prob_upars`, and
`log_lik_i_upars`. All of these functions should take `...`. as an
argument in addition to those specified for each function.

## Methods (by class)

- `loo_moment_match(default)`: A default method that takes as arguments
  a user-specified model object `x`, a `loo` object and user-specified
  functions `post_draws`, `log_lik_i`, `unconstrain_pars`,
  `log_prob_upars`, and `log_lik_i_upars`.

## References

Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2021).
Implicitly adaptive importance sampling. *Statistics and Computing*, 31,
16. doi:10.1007/s11222-020-09982-2. arXiv preprint arXiv:1906.08850.

## See also

[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md),
[`loo_moment_match_split()`](https://mc-stan.org/loo/dev/reference/loo_moment_match_split.md)

## Examples

``` r
# See the vignette for loo_moment_match()
```
