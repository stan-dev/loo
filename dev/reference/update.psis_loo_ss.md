# Update `psis_loo_ss` objects

Update `psis_loo_ss` objects

## Usage

``` r
# S3 method for class 'psis_loo_ss'
update(
  object,
  ...,
  data = NULL,
  draws = NULL,
  observations = NULL,
  r_eff = 1,
  cores = getOption("mc.cores", 1),
  loo_approximation = NULL,
  loo_approximation_draws = NULL,
  llgrad = NULL,
  llhess = NULL
)
```

## Arguments

- object:

  A `psis_loo_ss` object to update.

- ...:

  Currently not used.

- data, draws:

  See
  [`loo_subsample.function()`](https://mc-stan.org/loo/dev/reference/loo_subsample.md).

- observations:

  The subsample observations to use. The argument can take four (4)
  types of arguments:

  - `NULL` to use all observations. The algorithm then just uses
    standard [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) or
    [`loo_approximate_posterior()`](https://mc-stan.org/loo/dev/reference/loo_approximate_posterior.md).

  - A single integer to specify the number of observations to be
    subsampled.

  - A vector of integers to provide the indices used to subset the data.
    *These observations need to be subsampled with the same scheme as
    given by the `estimator` argument*.

  - A `psis_loo_ss` object to use the same observations that were used
    in a previous call to
    [`loo_subsample()`](https://mc-stan.org/loo/dev/reference/loo_subsample.md).

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

- loo_approximation:

  What type of approximation of the loo_i's should be used? The default
  is `"plpd"` (the log predictive density using the posterior
  expectation). There are six different methods implemented to
  approximate loo_i's (see the references for more details):

  - `"plpd"`: uses the lpd based on point estimates (i.e.,
    \\p(y_i\|\hat{\theta})\\).

  - `"lpd"`: uses the lpds (i,e., \\p(y_i\|y)\\).

  - `"tis"`: uses truncated importance sampling to approximate PSIS-LOO.

  - `"waic"`: uses waic (i.e., \\p(y_i\|y) - p\_{waic}\\).

  - `"waic_grad_marginal"`: uses waic approximation using first order
    delta method and posterior marginal variances to approximate
    \\p\_{waic}\\ (ie. \\p(y_i\|\hat{\theta})\\-p_waic_grad_marginal).
    Requires gradient of likelihood function.

  - `"waic_grad"`: uses waic approximation using first order delta
    method and posterior covariance to approximate \\p\_{waic}\\ (ie.
    \\p(y_i\|\hat{\theta})\\-p_waic_grad). Requires gradient of
    likelihood function.

  - `"waic_hess"`: uses waic approximation using second order delta
    method and posterior covariance to approximate \\p\_{waic}\\ (ie.
    \\p(y_i\|\hat{\theta})\\-p_waic_grad). Requires gradient and Hessian
    of likelihood function.

  As point estimates of \\\hat{\theta}\\, the posterior expectations of
  the parameters are used.

- loo_approximation_draws:

  The number of posterior draws used when integrating over the
  posterior. This is used if `loo_approximation` is set to `"lpd"`,
  `"waic"`, or `"tis"`.

- llgrad:

  The gradient of the log-likelihood. This is only used when
  `loo_approximation` is `"waic_grad"`, `"waic_grad_marginal"`, or
  `"waic_hess"`. The default is `NULL`.

- llhess:

  The Hessian of the log-likelihood. This is only used with
  `loo_approximation = "waic_hess"`. The default is `NULL`.

## Value

A `psis_loo_ss` object.

## Details

If `observations` is updated then if a vector of indices or a
`psis_loo_ss` object is supplied the updated object will have exactly
the observations indicated by the vector or `psis_loo_ss` object. If a
single integer is supplied, new observations will be sampled to reach
the supplied sample size.
