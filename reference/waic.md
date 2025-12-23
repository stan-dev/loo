# Widely applicable information criterion (WAIC)

The `waic()` methods can be used to compute WAIC from the pointwise
log-likelihood. However, we recommend LOO-CV using PSIS (as implemented
by the [`loo()`](https://mc-stan.org/loo/reference/loo.md) function)
because PSIS provides useful diagnostics as well as effective sample
size and Monte Carlo estimates.

## Usage

``` r
waic(x, ...)

# S3 method for class 'array'
waic(x, ...)

# S3 method for class 'matrix'
waic(x, ...)

# S3 method for class '`function`'
waic(x, ..., data = NULL, draws = NULL)

is.waic(x)
```

## Arguments

- x:

  A log-likelihood array, matrix, or function. The **Methods (by
  class)** section, below, has detailed descriptions of how to specify
  the inputs for each method.

- draws, data, ...:

  For the function method only. See the **Methods (by class)** section
  below for details on these arguments.

## Value

A named list (of class `c("waic", "loo")`) with components:

- `estimates`:

  A matrix with two columns (`"Estimate"`, `"SE"`) and three rows
  (`"elpd_waic"`, `"p_waic"`, `"waic"`). This contains point estimates
  and standard errors of the expected log pointwise predictive density
  (`elpd_waic`), the effective number of parameters (`p_waic`) and the
  information criterion `waic` (which is just `-2 * elpd_waic`, i.e.,
  converted to deviance scale).

- `pointwise`:

  A matrix with three columns (and number of rows equal to the number of
  observations) containing the pointwise contributions of each of the
  above measures (`elpd_waic`, `p_waic`, `waic`).

## Methods (by class)

- `waic(array)`: An \\I\\ by \\C\\ by \\N\\ array, where \\I\\ is the
  number of MCMC iterations per chain, \\C\\ is the number of chains,
  and \\N\\ is the number of data points.

- `waic(matrix)`: An \\S\\ by \\N\\ matrix, where \\S\\ is the size of
  the posterior sample (with all chains merged) and \\N\\ is the number
  of data points.

- `` waic(`function`) ``: A function `f()` that takes arguments `data_i`
  and `draws` and returns a vector containing the log-likelihood for a
  single observation `i` evaluated at each posterior draw. The function
  should be written such that, for each observation `i` in `1:N`,
  evaluating

      f(data_i = data[i,, drop=FALSE], draws = draws)

  results in a vector of length `S` (size of posterior sample). The
  log-likelihood function can also have additional arguments but
  `data_i` and `draws` are required.

  If using the function method then the arguments `data` and `draws`
  must also be specified in the call to
  [`loo()`](https://mc-stan.org/loo/reference/loo.md):

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

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely application information criterion in singular learning
theory. *Journal of Machine Learning Research* **11**, 3571-3594.

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

- The **loo** package [vignettes](https://mc-stan.org/loo/articles/) and
  Vehtari, Gelman, and Gabry (2017) and Vehtari, Simpson, Gelman, Yao,
  and Gabry (2024) for more details on why we prefer
  [`loo()`](https://mc-stan.org/loo/reference/loo.md) to `waic()`.

- [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.md)
  for comparing models on approximate LOO-CV or WAIC.

## Examples

``` r
### Array and matrix methods
LLarr <- example_loglik_array()
dim(LLarr)
#> [1] 500   2  32

LLmat <- example_loglik_matrix()
dim(LLmat)
#> [1] 1000   32

waic_arr <- waic(LLarr)
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
waic_mat <- waic(LLmat)
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
identical(waic_arr, waic_mat)
#> [1] TRUE


# \dontrun{
log_lik1 <- extract_log_lik(stanfit1)
#> Error: object 'stanfit1' not found
log_lik2 <- extract_log_lik(stanfit2)
#> Error: object 'stanfit2' not found
(waic1 <- waic(log_lik1))
#> Error: object 'log_lik1' not found
(waic2 <- waic(log_lik2))
#> Error: object 'log_lik2' not found
print(compare(waic1, waic2), digits = 2)
#> Warning: 'compare' is deprecated.
#> Use 'loo_compare' instead.
#> See help("Deprecated")
#> Error: object 'waic1' not found
# }
```
