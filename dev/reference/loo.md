# Efficient approximate leave-one-out cross-validation (LOO)

The `loo()` methods for arrays, matrices, and functions compute PSIS-LOO
CV, efficient approximate leave-one-out (LOO) cross-validation for
Bayesian models using Pareto smoothed importance sampling
([PSIS](https://mc-stan.org/loo/dev/reference/psis.md)). This is an
implementation of the methods described in Vehtari, Gelman, and Gabry
(2017) and Vehtari, Simpson, Gelman, Yao, and Gabry (2024).

The `loo_i()` function enables testing log-likelihood functions for use
with the `loo.function()` method.

## Usage

``` r
loo(x, ...)

# S3 method for class 'array'
loo(
  x,
  ...,
  r_eff = 1,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1),
  is_method = c("psis", "tis", "sis")
)

# S3 method for class 'matrix'
loo(
  x,
  ...,
  r_eff = 1,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1),
  is_method = c("psis", "tis", "sis")
)

# S3 method for class '`function`'
loo(
  x,
  ...,
  data = NULL,
  draws = NULL,
  r_eff = 1,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1),
  is_method = c("psis", "tis", "sis")
)

loo_i(i, llfun, ..., data = NULL, draws = NULL, r_eff = 1, is_method = "psis")

is.loo(x)

is.psis_loo(x)
```

## Arguments

- x:

  A log-likelihood array, matrix, or function. The **Methods (by
  class)** section, below, has detailed descriptions of how to specify
  the inputs for each method.

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

  Should the `psis` object created internally by `loo()` be saved in the
  returned object? The `loo()` function calls
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

- data, draws, ...:

  For the `loo.function()` method and the `loo_i()` function, these are
  the data, posterior draws, and other arguments to pass to the
  log-likelihood function. See the **Methods (by class)** section below
  for details on how to specify these arguments.

- i:

  For `loo_i()`, an integer in `1:N`.

- llfun:

  For `loo_i()`, the same as `x` for the `loo.function()` method. A
  log-likelihood function as described in the **Methods (by class)**
  section.

## Value

The `loo()` methods return a named list with class
`c("psis_loo", "loo")` and components:

- `estimates`:

  A matrix with two columns (`Estimate`, `SE`) and three rows
  (`elpd_loo`, `p_loo`, `looic`). This contains point estimates and
  standard errors of the expected log pointwise predictive density
  ([`elpd_loo`](https://mc-stan.org/loo/dev/reference/loo-glossary.md)),
  the effective number of parameters
  ([`p_loo`](https://mc-stan.org/loo/dev/reference/loo-glossary.md)) and
  the LOO information criterion `looic` (which is just `-2 * elpd_loo`,
  i.e., converted to deviance scale).

- `pointwise`:

  A matrix with five columns (and number of rows equal to the number of
  observations) containing the pointwise contributions of the measures
  (`elpd_loo`, `mcse_elpd_loo`, `p_loo`, `looic`, `influence_pareto_k`).
  in addition to the three measures in `estimates`, we also report
  pointwise values of the Monte Carlo standard error of
  [`elpd_loo`](https://mc-stan.org/loo/dev/reference/loo-glossary.md)
  ([`mcse_elpd_loo`](https://mc-stan.org/loo/dev/reference/loo-glossary.md)),
  and statistics describing the influence of each observation on the
  posterior distribution (`influence_pareto_k`). These are the estimates
  of the shape parameter \\k\\ of the generalized Pareto fit to the
  importance ratios for each leave-one-out distribution (see the
  [pareto-k-diagnostic](https://mc-stan.org/loo/dev/reference/pareto-k-diagnostic.md)
  page for details).

- `diagnostics`:

  A named list containing two vectors:

  - `pareto_k`: Importance sampling reliability diagnostics. By default,
    these are equal to the `influence_pareto_k` in `pointwise`. Some
    algorithms can improve importance sampling reliability and modify
    these diagnostics. See the
    [pareto-k-diagnostic](https://mc-stan.org/loo/dev/reference/pareto-k-diagnostic.md)
    page for details.

  - `n_eff`: PSIS effective sample size estimates.

- `psis_object`:

  This component will be `NULL` unless the `save_psis` argument is set
  to `TRUE` when calling `loo()`. In that case `psis_object` will be the
  object of class `"psis"` that is created when the `loo()` function
  calls [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md)
  internally to do the PSIS procedure.

The `loo_i()` function returns a named list with components `pointwise`
and `diagnostics`. These components have the same structure as the
`pointwise` and `diagnostics` components of the object returned by
`loo()` except they contain results for only a single observation.

## Details

The `loo()` function is an S3 generic and methods are provided for 3-D
pointwise log-likelihood arrays, pointwise log-likelihood matrices, and
log-likelihood functions. The array and matrix methods are the most
convenient, but for models fit to very large datasets the
`loo.function()` method is more memory efficient and may be preferable.

## Methods (by class)

- `loo(array)`: An \\I\\ by \\C\\ by \\N\\ array, where \\I\\ is the
  number of MCMC iterations per chain, \\C\\ is the number of chains,
  and \\N\\ is the number of data points.

- `loo(matrix)`: An \\S\\ by \\N\\ matrix, where \\S\\ is the size of
  the posterior sample (with all chains merged) and \\N\\ is the number
  of data points.

- `` loo(`function`) ``: A function `f()` that takes arguments `data_i`
  and `draws` and returns a vector containing the log-likelihood for a
  single observation `i` evaluated at each posterior draw. The function
  should be written such that, for each observation `i` in `1:N`,
  evaluating

      f(data_i = data[i,, drop=FALSE], draws = draws)

  results in a vector of length `S` (size of posterior sample). The
  log-likelihood function can also have additional arguments but
  `data_i` and `draws` are required.

  If using the function method then the arguments `data` and `draws`
  must also be specified in the call to `loo()`:

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

## Defining `loo()` methods in a package

Package developers can define `loo()` methods for fitted models objects.
See the example `loo.stanfit()` method in the **Examples** section below
for an example of defining a method that calls `loo.array()`. The
`loo.stanreg()` method in the **rstanarm** package is an example of
defining a method that calls `loo.function()`.

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

- The **loo** package
  [vignettes](https://mc-stan.org/loo/articles/index.html) for
  demonstrations.

- The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html)
  on the **loo** website for answers to frequently asked questions.

- [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md) for the
  underlying Pareto Smoothed Importance Sampling (PSIS) procedure used
  in the LOO-CV approximation.

- [pareto-k-diagnostic](https://mc-stan.org/loo/dev/reference/pareto-k-diagnostic.md)
  for convenience functions for looking at diagnostics.

- [`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md)
  for model comparison.

## Examples

``` r
### Array and matrix methods (using example objects included with loo package)
# Array method
LLarr <- example_loglik_array()
rel_n_eff <- relative_eff(exp(LLarr))
loo(LLarr, r_eff = rel_n_eff, cores = 2)
#> 
#> Computed from 1000 by 32 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -83.6 4.3
#> p_loo         3.3 1.2
#> looic       167.2 8.6
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.0]).
#> 
#> All Pareto k estimates are good (k < 0.67).
#> See help('pareto-k-diagnostic') for details.

# Matrix method
LLmat <- example_loglik_matrix()
rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = 500))
loo(LLmat, r_eff = rel_n_eff, cores = 2)
#> 
#> Computed from 1000 by 32 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -83.6 4.3
#> p_loo         3.3 1.2
#> looic       167.2 8.6
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.0]).
#> 
#> All Pareto k estimates are good (k < 0.67).
#> See help('pareto-k-diagnostic') for details.


### Using log-likelihood function instead of array or matrix
set.seed(124)

# Simulate data and draw from posterior
N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
p <- rbeta(1, a0, b0)
y <- rbinom(N, size = K, prob = p)
a <- a0 + sum(y); b <- b0 + N * K - sum(y)
fake_posterior <- as.matrix(rbeta(S, a, b))
dim(fake_posterior) # S x 1
#> [1] 100   1
fake_data <- data.frame(y,K)
dim(fake_data) # N x 2
#> [1] 50  2

llfun <- function(data_i, draws) {
  # each time called internally within loo the arguments will be equal to:
  # data_i: ith row of fake_data (fake_data[i,, drop=FALSE])
  # draws: entire fake_posterior matrix
  dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
}

# Use the loo_i function to check that llfun works on a single observation
# before running on all obs. For example, using the 3rd obs in the data:
loo_3 <- loo_i(i = 3, llfun = llfun, data = fake_data, draws = fake_posterior)
print(loo_3$pointwise[, "elpd_loo"])
#>  elpd_loo 
#> -1.267103 

# Use loo.function method (default r_eff=1 is used as this posterior not obtained via MCMC)
loo_with_fn <- loo(llfun, draws = fake_posterior, data = fake_data)

# If we look at the elpd_loo contribution from the 3rd obs it should be the
# same as what we got above with the loo_i function and i=3:
print(loo_with_fn$pointwise[3, "elpd_loo"])
#>  elpd_loo 
#> -1.267103 
print(loo_3$pointwise[, "elpd_loo"])
#>  elpd_loo 
#> -1.267103 

# Check that the loo.matrix method gives same answer as loo.function method
log_lik_matrix <- sapply(1:N, function(i) {
  llfun(data_i = fake_data[i,, drop=FALSE], draws = fake_posterior)
})
loo_with_mat <- loo(log_lik_matrix)
all.equal(loo_with_mat$estimates, loo_with_fn$estimates) # should be TRUE!
#> [1] TRUE


# \dontrun{
### For package developers: defining loo methods

# An example of a possible loo method for 'stanfit' objects (rstan package).
# A similar method is included in the rstan package.
# In order for users to be able to call loo(stanfit) instead of
# loo.stanfit(stanfit) the NAMESPACE needs to be handled appropriately
# (roxygen2 and devtools packages are good for that).
#
loo.stanfit <-
 function(x,
         pars = "log_lik",
         ...,
         save_psis = FALSE,
         cores = getOption("mc.cores", 1)) {
  stopifnot(length(pars) == 1L)
  LLarray <- loo::extract_log_lik(stanfit = x,
                                  parameter_name = pars,
                                  merge_chains = FALSE)
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  loo::loo.array(LLarray,
                 r_eff = r_eff,
                 cores = cores,
                 save_psis = save_psis)
}
# }

```
