# Model comparison

Compare fitted models based on
[ELPD](https://mc-stan.org/loo/reference/loo-glossary.md).

By default the print method shows only the most important information.
Use `print(..., simplify=FALSE)` to print a more detailed summary.

## Usage

``` r
loo_compare(x, ...)

# Default S3 method
loo_compare(x, ...)

# S3 method for class 'compare.loo'
print(x, ..., digits = 1, simplify = TRUE)

# S3 method for class 'compare.loo_ss'
print(x, ..., digits = 1, simplify = TRUE)
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

- simplify:

  For the print method only, should only the essential columns of the
  summary matrix be printed? The entire matrix is always returned, but
  by default only the most important columns are printed.

## Value

A matrix with class `"compare.loo"` that has its own print method. See
the **Details** section.

## Details

When comparing two fitted models, we can estimate the difference in
their expected predictive accuracy by the difference in
[`elpd_loo`](https://mc-stan.org/loo/reference/loo-glossary.md) or
`elpd_waic` (or multiplied by \\-2\\, if desired, to be on the deviance
scale).

When using `loo_compare()`, the returned matrix will have one row per
model and several columns of estimates. The values in the
[`elpd_diff`](https://mc-stan.org/loo/reference/loo-glossary.md) and
[`se_diff`](https://mc-stan.org/loo/reference/loo-glossary.md) columns
of the returned matrix are computed by making pairwise comparisons
between each model and the model with the largest ELPD (the model in the
first row). For this reason the `elpd_diff` column will always have the
value `0` in the first row (i.e., the difference between the preferred
model and itself) and negative values in subsequent rows for the
remaining models.

To compute the standard error of the difference in
[ELPD](https://mc-stan.org/loo/reference/loo-glossary.md) — which should
not be expected to equal the difference of the standard errors — we use
a paired estimate to take advantage of the fact that the same set of
\\N\\ data points was used to fit both models. These calculations should
be most useful when \\N\\ is large, because then non-normality of the
distribution is not such an issue when estimating the uncertainty in
these sums. These standard errors, for all their flaws, should give a
better sense of uncertainty than what is obtained using the current
standard approach of comparing differences of deviances to a Chi-squared
distribution, a practice derived for Gaussian linear models or
asymptotically, and which only applies to nested models in any case.
Sivula et al. (2022) discuss the conditions when the normal
approximation used for SE and `se_diff` is good.

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
#>        elpd_diff se_diff
#> model3   0.00      0.00 
#> model2 -32.00      0.00 
#> model1 -64.00      0.00 

# show more details with simplify=FALSE
# (will be the same for all models in this artificial example)
print(comp, simplify = FALSE, digits = 3)
#>        elpd_diff se_diff elpd_loo se_elpd_loo p_loo   se_p_loo looic   se_looic
#> model3   0.000     0.000 -19.589    4.284       3.329   1.152   39.178   8.568 
#> model2 -32.000     0.000 -51.589    4.284       3.329   1.152  103.178   8.568 
#> model1 -64.000     0.000 -83.589    4.284       3.329   1.152  167.178   8.568 

# can use a list of objects with custom names
# will use apple, banana, and cherry, as the names in the output
loo_compare(list("apple" = loo1, "banana" = loo2, "cherry" = loo3))
#>        elpd_diff se_diff
#> cherry   0.0       0.0  
#> banana -32.0       0.0  
#> apple  -64.0       0.0  

# \dontrun{
# works for waic (and kfold) too
loo_compare(waic(LL), waic(LL - 10))
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#> Warning: 
#> 3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#>        elpd_diff se_diff
#> model1    0.0       0.0 
#> model2 -320.0       0.0 
# }
```
