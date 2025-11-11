# Efficient approximate leave-one-out cross-validation (LOO) using subsampling, so that less costly and more approximate computation is made for all LOO-fold, and more costly and accurate computations are made only for m\<N LOO-folds.

Efficient approximate leave-one-out cross-validation (LOO) using
subsampling, so that less costly and more approximate computation is
made for all LOO-fold, and more costly and accurate computations are
made only for m\<N LOO-folds.

## Usage

``` r
loo_subsample(x, ...)

# S3 method for class '`function`'
loo_subsample(
  x,
  ...,
  data = NULL,
  draws = NULL,
  observations = 400,
  log_p = NULL,
  log_g = NULL,
  r_eff = 1,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1),
  loo_approximation = "plpd",
  loo_approximation_draws = NULL,
  estimator = "diff_srs",
  llgrad = NULL,
  llhess = NULL
)
```

## Arguments

- x:

  A function. The **Methods (by class)** section, below, has detailed
  descriptions of how to specify the inputs.

- data, draws, ...:

  For `loo_subsample.function()`, these are the data, posterior draws,
  and other arguments to pass to the log-likelihood function. Note that
  for some `loo_approximation`s, the draws will be replaced by the
  posteriors summary statistics to compute loo approximations. See
  argument `loo_approximation` for details.

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
    in a previous call to `loo_subsample()`.

- log_p, log_g:

  Should be supplied only if approximate posterior draws are used. The
  default (`NULL`) indicates draws are from "true" posterior (i.e. using
  MCMC). If not `NULL` then they should be specified as described in
  [`loo_approximate_posterior()`](https://mc-stan.org/loo/dev/reference/loo_approximate_posterior.md).

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

  Should the `"psis"` object created internally by `loo_subsample()` be
  saved in the returned object? See
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) for details.

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

- estimator:

  How should `elpd_loo`, `p_loo` and `looic` be estimated? The default
  is `"diff_srs"`.

  - `"diff_srs"`: uses the difference estimator with simple random
    sampling without replacement (srs). `p_loo` is estimated using
    standard srs. (Magnusson et al., 2020)

  - `"hh"`: uses the Hansen-Hurwitz estimator with sampling with
    replacement proportional to size, where `abs` of loo_approximation
    is used as size. (Magnusson et al., 2019)

  - `"srs"`: uses simple random sampling and ordinary estimation.

- llgrad:

  The gradient of the log-likelihood. This is only used when
  `loo_approximation` is `"waic_grad"`, `"waic_grad_marginal"`, or
  `"waic_hess"`. The default is `NULL`.

- llhess:

  The Hessian of the log-likelihood. This is only used with
  `loo_approximation = "waic_hess"`. The default is `NULL`.

## Value

`loo_subsample()` returns a named list with class
`c("psis_loo_ss", "psis_loo", "loo")`. This has the same structure as
objects returned by
[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) but with the
additional slot:

- `loo_subsampling`: A list with two vectors, `log_p` and `log_g`, of
  the same length containing the posterior density and the approximation
  density for the individual draws.

## Details

The `loo_subsample()` function is an S3 generic and a methods is
currently provided for log-likelihood functions. The implementation
works for both MCMC and for posterior approximations where it is
possible to compute the log density for the approximation.

## Methods (by class)

- `` loo_subsample(`function`) ``: A function `f()` that takes arguments
  `data_i` and `draws` and returns a vector containing the
  log-likelihood for a single observation `i` evaluated at each
  posterior draw. The function should be written such that, for each
  observation `i` in `1:N`, evaluating

      f(data_i = data[i,, drop=FALSE], draws = draws)

  results in a vector of length `S` (size of posterior sample). The
  log-likelihood function can also have additional arguments but
  `data_i` and `draws` are required.

  If using the function method then the arguments `data` and `draws`
  must also be specified in the call to
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md):

  - `data`: A data frame or matrix containing the data (e.g. observed
    outcome and predictors) needed to compute the pointwise
    log-likelihood. For each observation `i`, the `i`th row of `data`
    will be passed to the `data_i` argument of the log-likelihood
    function.

  - `draws`: An object containing the posterior draws for any parameters
    needed to compute the pointwise log-likelihood. Unlike `data`, which
    is indexed by observation, for each observation the entire object
    `draws` will be passed to the `draws` argument of the log-likelihood
    function.

  - The `...` can be used if your log-likelihood function takes
    additional arguments. These arguments are used like the `draws`
    argument in that they are recycled for each observation.

## References

Magnusson, M., Riis Andersen, M., Jonasson, J. and Vehtari, A. (2019).
Leave-One-Out Cross-Validation for Large Data. In *Thirty-sixth
International Conference on Machine Learning*, PMLR 97:4244-4253.

Magnusson, M., Riis Andersen, M., Jonasson, J. and Vehtari, A. (2020).
Leave-One-Out Cross-Validation for Model Comparison in Large Data. In
*Proceedings of the 23rd International Conference on Artificial
Intelligence and Statistics (AISTATS)*, PMLR 108:341-351.

## See also

[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md),
[`psis()`](https://mc-stan.org/loo/dev/reference/psis.md),
[`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md)
