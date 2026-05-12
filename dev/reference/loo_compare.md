# Model comparison

Compare fitted models based on
[ELPD](https://mc-stan.org/loo/dev/reference/loo-glossary.md).

## Usage

``` r
loo_compare(x, ...)

# Default S3 method
loo_compare(x, ...)

# S3 method for class 'compare.loo'
print(x, ..., digits = 1, p_worse = TRUE)

# S3 method for class 'compare.loo_ss'
print(x, ..., digits = 1)
```

## Arguments

- x:

  An object of class `"loo"` or a list of such objects. If a list is
  used then the list names will be used as the model names in the
  output. See **Examples**.

- ...:

  Additional objects of class `"loo"`, if not passed in as a single
  list.

- digits:

  For the print method only, the number of digits to use when printing.

- p_worse:

  For the print method only, should we include the normal approximation
  based probability of each model having worse performance than the best
  model? The default is `TRUE`.

## Value

A data frame with class `"compare.loo"` that has its own print method.
See the **Details** and **Examples** sections.

## Details

When comparing two fitted models, we can estimate the difference in
their expected predictive accuracy by the difference in
[`elpd_loo`](https://mc-stan.org/loo/dev/reference/loo-glossary.md) or
`elpd_waic` (or multiplied by \\-2\\, if desired, to be on the deviance
scale).

### `elpd_diff` and `se_diff`

When using `loo_compare()`, the returned data frame will have one row
per model and several columns of estimates. The values of
[`elpd_diff`](https://mc-stan.org/loo/dev/reference/loo-glossary.md) and
[`se_diff`](https://mc-stan.org/loo/dev/reference/loo-glossary.md) are
computed by making pairwise comparisons between each model and the model
with the largest ELPD (the model listed first). Therefore, the first
`elpd_diff` value will always be `0` (i.e., the difference between the
preferred model and itself) and the rest of the values will be negative.

To compute the standard error of the difference in
[ELPD](https://mc-stan.org/loo/dev/reference/loo-glossary.md) — which
should not be expected to equal the difference of the standard errors —
we use a paired estimate to take advantage of the fact that the same set
of \\N\\ data points was used to fit both models. These calculations
should be most useful when \\N\\ is large, because then non-normality of
the distribution is not such an issue when estimating the uncertainty in
these sums. These standard errors, for all their flaws, should give a
better sense of uncertainty than what is obtained using the current
standard approach of comparing differences of deviances to a Chi-squared
distribution, a practice derived for Gaussian linear models or
asymptotically, and which only applies to nested models in any case.

### `p_worse`, `diag_diff`, and `diag_elpd`

The values in the `p_worse` column show the probability of each model
having worse ELPD than the best model. These probabilities are computed
with a normal approximation using the values from `elpd_diff` and
`se_diff`. Sivula et al. (2025) present the conditions when the normal
approximation used for SE and `se_diff` is good, and the column
`diag_diff` contains possible diagnostic messages:

- `N < 100` (small data)

- `|elpd_diff| < 4` (models make similar predictions)

If either of these diagnostic messages is shown, the error distribution
is skewed or thick tailed and the normal approximation based on
`elpd_diff` and `se_diff` is not well calibrated. In that case, the
probabilities `p_worse` are likely to be too large. However, `elpd_diff`
and `se_diff` will still be indicative of the differences and
uncertainties (for example, if `|elpd_diff|` is many times larger than
`se_diff` the difference is quite certain). In addition, if the model is
not well specificed and there are outliers, the error distribution can
also be skewed or thick tailed and the normal approximation is not well
calibrated. Possible model misspecification and outliers can be
diagnosed with usual predictive checking methods.

The column `diag_elpd` shows the PSIS-LOO Pareto k diagnostic for the
pointwise ELPD computations for each model. If `K k_psis > 0.7` is
shown, where `K` is the number of high Pareto k values in the PSIS
computation, then there may be significant bias in `elpd_diff` favoring
models with a large number of high Pareto k values.

### Warnings for many model comparisons

If more than \\11\\ models are compared, we internally recompute the
model differences using the median model by ELPD as the baseline model.
We then estimate whether the differences in predictive performance are
potentially due to chance as described by McLatchie and Vehtari (2023).
This will flag a warning if it is deemed that there is a risk of
over-fitting due to the selection process. In that case users are
recommended to avoid model selection based on LOO-CV, and instead to
favor model averaging/stacking or projection predictive inference.

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

Sivula, T, Magnusson, M., Matamoros A. A., and Vehtari, A. (2025).
Uncertainty in Bayesian leave-one-out cross-validation based model
comparison. *Bayesian Analysis*.
[doi:10.1214/25-BA1569](https://doi.org/10.1214/25-BA1569)

McLatchie, Y., and Vehtari, A. (2024). Efficient estimation and
correction of selection-induced bias with order statistics. *Statistics
and Computing*. 34(132).
[doi:10.1007/s11222-024-10442-4](https://doi.org/10.1007/s11222-024-10442-4)

## See also

- The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html)
  on the **loo** website for answers to frequently asked questions.

## Examples

``` r
# very artificial example, just for demonstration!
LL <- example_loglik_array()
loo1 <- loo(LL)     # should be worst model when compared
loo2 <- loo(LL + 1) # should be second best model when compared
loo3 <- loo(LL + 2) # should be best model when compared

comp <- loo_compare(loo1, loo2, loo3)
print(comp, digits = 2)
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model3      0.00    0.00      NA                    
#>  model2    -32.00    0.00    1.00   N < 100          
#>  model1    -64.00    0.00    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.

# can use a list of objects with custom names
# the names will be used in the output
loo_compare(list("apple" = loo1, "banana" = loo2, "cherry" = loo3))
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  cherry       0.0     0.0      NA                    
#>  banana     -32.0     0.0    1.00   N < 100          
#>   apple     -64.0     0.0    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.

# \dontrun{
# works for waic (and kfold) too
loo_compare(waic(LL), waic(LL - 10))
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model1       0.0     0.0      NA                    
#>  model2    -320.0     0.0    1.00   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
