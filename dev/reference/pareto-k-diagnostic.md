# Diagnostics for Pareto smoothed importance sampling (PSIS)

Print a diagnostic table summarizing the estimated Pareto shape
parameters and PSIS effective sample sizes, find the indexes of
observations for which the estimated Pareto shape parameter \\k\\ is
larger than some `threshold` value, or plot observation indexes vs.
diagnostic estimates. The **Details** section below provides a brief
overview of the diagnostics, but we recommend consulting Vehtari,
Gelman, and Gabry (2017) and Vehtari, Simpson, Gelman, Yao, and Gabry
(2024) for full details.

## Usage

``` r
pareto_k_table(x)

pareto_k_ids(x, threshold = NULL)

pareto_k_values(x)

pareto_k_influence_values(x)

psis_n_eff_values(x)

mcse_loo(x, threshold = NULL)

# S3 method for class 'psis_loo'
plot(
  x,
  diagnostic = c("k", "ESS", "n_eff"),
  ...,
  label_points = FALSE,
  main = "PSIS diagnostic plot"
)

# S3 method for class 'psis'
plot(
  x,
  diagnostic = c("k", "ESS", "n_eff"),
  ...,
  label_points = FALSE,
  main = "PSIS diagnostic plot"
)
```

## Arguments

- x:

  An object created by
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) or
  [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md).

- threshold:

  For `pareto_k_ids()`, `threshold` is the minimum \\k\\ value to flag
  (default is a sample size `S` dependend threshold `1 - 1 / log10(S)`).
  For `mcse_loo()`, if any \\k\\ estimates are greater than `threshold`
  the MCSE estimate is returned as `NA` See **Details** for the
  motivation behind these defaults.

- diagnostic:

  For the `plot` method, which diagnostic should be plotted? The options
  are `"k"` for Pareto \\k\\ estimates (the default), or `"ESS"` or
  `"n_eff"` for PSIS effective sample size estimates.

- label_points, ...:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method, if `label_points` is `TRUE` the observation numbers
  corresponding to any values of \\k\\ greater than the diagnostic
  threshold will be displayed in the plot. Any arguments specified in
  `...` will be passed to
  [`graphics::text()`](https://rdrr.io/r/graphics/text.html) and can be
  used to control the appearance of the labels.

- main:

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method, a title for the plot.

## Value

`pareto_k_table()` returns an object of class `"pareto_k_table"`, which
is a matrix with columns `"Count"`, `"Proportion"`, and `"Min. n_eff"`,
and has its own print method.

`pareto_k_ids()` returns an integer vector indicating which observations
have Pareto \\k\\ estimates above `threshold`.

`pareto_k_values()` returns a vector of the estimated Pareto \\k\\
parameters. These represent the reliability of sampling.

`pareto_k_influence_values()` returns a vector of the estimated Pareto
\\k\\ parameters. These represent influence of the observations on the
model posterior distribution.

`psis_n_eff_values()` returns a vector of the estimated PSIS effective
sample sizes.

`mcse_loo()` returns the Monte Carlo standard error (MCSE) estimate for
PSIS-LOO. MCSE will be NA if any Pareto \\k\\ values are above
`threshold`.

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method is
called for its side effect and does not return anything. If `x` is the
result of a call to
[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) or
[`psis()`](https://mc-stan.org/loo/dev/reference/psis.md) then
`plot(x, diagnostic)` produces a plot of the estimates of the Pareto
shape parameters (`diagnostic = "k"`) or estimates of the PSIS effective
sample sizes (`diagnostic = "ESS"`).

## Details

The reliability and approximate convergence rate of the PSIS-based
estimates can be assessed using the estimates for the shape parameter
\\k\\ of the generalized Pareto distribution. The diagnostic threshold
for Pareto \\k\\ depends on sample size \\S\\ (sample size dependent
threshold was introduced by Vehtari et al. (2024), and before that fixed
thresholds of 0.5 and 0.7 were recommended). For simplicity, `loo`
package uses the nominal sample size \\S\\ when computing the sample
size specific threshold. This provides an optimistic threshold if the
effective sample size is less than 2200, but if MCMC-ESS \> S/2 the
difference is usually negligible. Thinning of MCMC draws can be used to
improve the ratio ESS/S.

- If \\k \< min(1 - 1 / log10(S), 0.7)\\, where \\S\\ is the sample
  size, the PSIS estimate and the corresponding Monte Carlo standard
  error estimate are reliable.

- If \\1 - 1 / log10(S) \<= k \< 0.7\\, the PSIS estimate and the
  corresponding Monte Carlo standard error estimate are not reliable,
  but increasing the (effective) sample size \\S\\ above 2200 may help
  (this will increase the sample size specific threshold
  \\(1-1/log10(2200)\>0.7\\ and then the bias specific threshold 0.7
  dominates).

- If \\0.7 \<= k \< 1\\, the PSIS estimate and the corresponding Monte
  Carlo standard error have large bias and are not reliable. Increasing
  the sample size may reduce the variability in \\k\\ estimate, which
  may result in lower \\k\\ estimate, too.

- If \\k \geq 1\\, the target distribution is estimated to have a
  non-finite mean. The PSIS estimate and the corresponding Monte Carlo
  standard error are not well defined. Increasing the sample size may
  reduce the variability in the \\k\\ estimate, which may also result in
  a lower \\k\\ estimate.

### What if the estimated tail shape parameter \\k\\ exceeds the diagnostic threshold?

Importance sampling is likely to work less well if the marginal
posterior \\p(\theta^s \| y)\\ and LOO posterior \\p(\theta^s \|
y\_{-i})\\ are very different, which is more likely to happen with a
non-robust model and highly influential observations. If the estimated
tail shape parameter \\k\\ exceeds the diagnostic threshold, the user
should be warned. (Note: If \\k\\ is greater than the diagnostic
threshold then WAIC is also likely to fail, but WAIC lacks as accurate
diagnostic.) When using PSIS in the context of approximate LOO-CV, we
recommend one of the following actions:

- With some additional computations, it is possible to transform the
  MCMC draws from the posterior distribution to obtain more reliable
  importance sampling estimates. This results in a smaller shape
  parameter \\k\\. See
  [`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
  and the vignette *Avoiding model refits in leave-one-out
  cross-validation with moment matching* for an example of this.

- Sampling from a leave-one-out mixture distribution (see the vignette
  *Mixture IS leave-one-out cross-validation for high-dimensional
  Bayesian models*), directly from \\p(\theta^s \| y\_{-i})\\ for the
  problematic observations \\i\\, or using \\K\\-fold cross-validation
  (see the vignette *Holdout validation and K-fold cross-validation of
  Stan programs with the loo package*) will generally be more stable.

- Using a model that is more robust to anomalous observations will
  generally make approximate LOO-CV more stable.

### Observation influence statistics

The estimated shape parameter \\k\\ for each observation can be used as
a measure of the observation's influence on posterior distribution of
the model. These can be obtained with `pareto_k_influence_values()`.

### Effective sample size and error estimates

In the case that we obtain the samples from the proposal distribution
via MCMC the **loo** package also computes estimates for the Monte Carlo
error and the effective sample size for importance sampling, which are
more accurate for PSIS than for IS and TIS (see Vehtari et al (2024) for
details). However, the PSIS effective sample size estimate will be
**over-optimistic when the estimate of \\k\\ is greater than**
\\min(1-1/log10(S), 0.7)\\, where \\S\\ is the sample size.

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

- [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md) for the
  implementation of the PSIS algorithm.

- The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html)
  on the **loo** website for answers to frequently asked questions.
