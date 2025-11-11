# Model comparison (deprecated, old version)

**This function is deprecated**. Please use the new
[`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md)
function instead.

## Usage

``` r
compare(..., x = list())
```

## Arguments

- ...:

  At least two objects returned by
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) (or
  [`waic()`](https://mc-stan.org/loo/dev/reference/waic.md)).

- x:

  A list of at least two objects returned by
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) (or
  [`waic()`](https://mc-stan.org/loo/dev/reference/waic.md)). This
  argument can be used as an alternative to specifying the objects in
  `...`.

## Value

A vector or matrix with class `'compare.loo'` that has its own print
method. If exactly two objects are provided in `...` or `x`, then the
difference in expected predictive accuracy and the standard error of the
difference are returned. If more than two objects are provided then a
matrix of summary information is returned (see **Details**).

## Details

When comparing two fitted models, we can estimate the difference in
their expected predictive accuracy by the difference in `elpd_loo` or
`elpd_waic` (or multiplied by -2, if desired, to be on the deviance
scale).

*When that difference, `elpd_diff`, is positive then the expected
predictive accuracy for the second model is higher. A negative
`elpd_diff` favors the first model.*

When using `compare()` with more than two models, the values in the
`elpd_diff` and `se_diff` columns of the returned matrix are computed by
making pairwise comparisons between each model and the model with the
best ELPD (i.e., the model in the first row). Although the `elpd_diff`
column is equal to the difference in `elpd_loo`, do not expect the
`se_diff` column to be equal to the the difference in `se_elpd_loo`.

To compute the standard error of the difference in ELPD we use a paired
estimate to take advantage of the fact that the same set of *N* data
points was used to fit both models. These calculations should be most
useful when *N* is large, because then non-normality of the distribution
is not such an issue when estimating the uncertainty in these sums.
These standard errors, for all their flaws, should give a better sense
of uncertainty than what is obtained using the current standard approach
of comparing differences of deviances to a Chi-squared distribution, a
practice derived for Gaussian linear models or asymptotically, and which
only applies to nested models in any case.

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

## Examples

``` r
# \dontrun{
loo1 <- loo(log_lik1)
#> Error: object 'log_lik1' not found
loo2 <- loo(log_lik2)
#> Error: object 'log_lik2' not found
print(compare(loo1, loo2), digits = 3)
#> Warning: 'compare' is deprecated.
#> Use 'loo_compare' instead.
#> See help("Deprecated")
#> Error: object 'loo1' not found
print(compare(x = list(loo1, loo2)))
#> Warning: 'compare' is deprecated.
#> Use 'loo_compare' instead.
#> See help("Deprecated")
#> Error: object 'loo1' not found

waic1 <- waic(log_lik1)
#> Error: object 'log_lik1' not found
waic2 <- waic(log_lik2)
#> Error: object 'log_lik2' not found
compare(waic1, waic2)
#> Warning: 'compare' is deprecated.
#> Use 'loo_compare' instead.
#> See help("Deprecated")
#> Error: object 'waic1' not found
# }
```
