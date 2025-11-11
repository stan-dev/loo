# Efficient approximate leave-one-out cross-validation (LOO) for posterior approximations

Efficient approximate leave-one-out cross-validation (LOO) for posterior
approximations

## Usage

``` r
loo_approximate_posterior(x, log_p, log_g, ...)

# S3 method for class 'array'
loo_approximate_posterior(
  x,
  log_p,
  log_g,
  ...,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1)
)

# S3 method for class 'matrix'
loo_approximate_posterior(
  x,
  log_p,
  log_g,
  ...,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1)
)

# S3 method for class '`function`'
loo_approximate_posterior(
  x,
  ...,
  data = NULL,
  draws = NULL,
  log_p = NULL,
  log_g = NULL,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1)
)
```

## Arguments

- x:

  A log-likelihood array, matrix, or function. The **Methods (by
  class)** section, below, has detailed descriptions of how to specify
  the inputs for each method.

- log_p:

  The log-posterior (target) evaluated at S samples from the proposal
  distribution (g). A vector of length S.

- log_g:

  The log-density (proposal) evaluated at S samples from the proposal
  distribution (g). A vector of length S.

- save_psis:

  Should the `"psis"` object created internally by
  `loo_approximate_posterior()` be saved in the returned object? See
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

- data, draws, ...:

  For the `loo_approximate_posterior.function()` method, these are the
  data, posterior draws, and other arguments to pass to the
  log-likelihood function. See the **Methods (by class)** section below
  for details on how to specify these arguments.

## Value

The `loo_approximate_posterior()` methods return a named list with class
`c("psis_loo_ap", "psis_loo", "loo")`. It has the same structure as the
objects returned by
[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) but with the
additional slot:

- `posterior_approximation`:

  A list with two vectors, `log_p` and `log_g` of the same length
  containing the posterior density and the approximation density for the
  individual draws.

## Details

The `loo_approximate_posterior()` function is an S3 generic and methods
are provided for 3-D pointwise log-likelihood arrays, pointwise
log-likelihood matrices, and log-likelihood functions. The
implementation works for posterior approximations where it is possible
to compute the log density for the posterior approximation.

## Methods (by class)

- `loo_approximate_posterior(array)`: An \\I\\ by \\C\\ by \\N\\ array,
  where \\I\\ is the number of MCMC iterations per chain, \\C\\ is the
  number of chains, and \\N\\ is the number of data points.

- `loo_approximate_posterior(matrix)`: An \\S\\ by \\N\\ matrix, where
  \\S\\ is the size of the posterior sample (with all chains merged) and
  \\N\\ is the number of data points.

- `` loo_approximate_posterior(`function`) ``: A function `f()` that
  takes arguments `data_i` and `draws` and returns a vector containing
  the log-likelihood for a single observation `i` evaluated at each
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
