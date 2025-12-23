# LOO package glossary

The pages provides definitions to key terms. Also see the [FAQ
page](https://mc-stan.org/loo/articles/online-only/faq.html) on the
**loo** website for answers to frequently asked questions.

Note: VGG2017 refers to Vehtari, Gelman, and Gabry (2017). See
**References**, below.

## ELPD and `elpd_loo`

The ELPD is the theoretical expected log pointwise predictive density
for a new dataset (Eq 1 in VGG2017), which can be estimated, e.g., using
cross-validation. `elpd_loo` is the Bayesian LOO estimate of the
expected log pointwise predictive density (Eq 4 in VGG2017) and is a sum
of N individual pointwise log predictive densities. Probability
densities can be smaller or larger than 1, and thus log predictive
densities can be negative or positive. For simplicity the ELPD acronym
is used also for expected log pointwise predictive probabilities for
discrete models. Probabilities are always equal or less than 1, and thus
log predictive probabilities are 0 or negative.

## Standard error of `elpd_loo`

As `elpd_loo` is defined as the sum of N independent components (Eq 4 in
VGG2017), we can compute the standard error by using the standard
deviation of the N components and multiplying by `sqrt(N)` (Eq 23 in
VGG2017). This standard error is a coarse description of our uncertainty
about the predictive performance for unknown future data. When N is
small or there is severe model misspecification, the current SE estimate
is overoptimistic and the actual SE can even be twice as large. Even for
moderate N, when the SE estimate is an accurate estimate for the scale,
it ignores the skewness. When making model comparisons, the SE of the
component-wise (pairwise) differences should be used instead (see the
`se_diff` section below and Eq 24 in VGG2017). Sivula et al. (2022)
discuss the conditions when the normal approximation used for SE and
`se_diff` is good.

## Monte Carlo SE of elpd_loo

The Monte Carlo standard error is the estimate for the computational
accuracy of MCMC and importance sampling used to compute `elpd_loo`.
Usually this is negligible compared to the standard describing the
uncertainty due to finite number of observations (Eq 23 in VGG2017).

## `p_loo` (effective number of parameters)

`p_loo` is the difference between `elpd_loo` and the non-cross-validated
log posterior predictive density. It describes how much more difficult
it is to predict future data than the observed data. Asymptotically
under certain regularity conditions, `p_loo` can be interpreted as the
*effective number of parameters*. In well behaving cases `p_loo < N` and
`p_loo < p`, where `p` is the total number of parameters in the model.
`p_loo > N` or `p_loo > p` indicates that the model has very weak
predictive capability and may indicate a severe model misspecification.
See below for more on interpreting `p_loo` when there are warnings about
high Pareto k diagnostic values.

## Pareto k estimates

The Pareto \\k\\ estimate is a diagnostic for Pareto smoothed importance
sampling (PSIS), which is used to compute components of `elpd_loo`. In
importance-sampling LOO the full posterior distribution is used as the
proposal distribution. The Pareto k diagnostic estimates how far an
individual leave-one-out distribution is from the full distribution. If
leaving out an observation changes the posterior too much then
importance sampling is not able to give a reliable estimate. Pareto
smoothing stabilizes importance sampling and guarantees a finite
variance estimate at the cost of some bias.

The diagnostic threshold for Pareto \\k\\ depends on sample size \\S\\
(sample size dependent threshold was introduced by Vehtari et al., 2024,
and before that fixed thresholds of 0.5 and 0.7 were recommended). For
simplicity, `loo` package uses the nominal sample size \\S\\ when
computing the sample size specific threshold. This provides an
optimistic threshold if the effective sample size is less than 2200, but
even then if ESS/S \> 1/2 the difference is usually negligible. Thinning
of MCMC draws can be used to improve the ratio ESS/S.

- If \\k \< \min(1 - 1 / \log\_{10}(S), 0.7)\\, where \\S\\ is the
  sample size, the PSIS estimate and the corresponding Monte Carlo
  standard error estimate are reliable.

- If \\1 - 1 / \log\_{10}(S) \<= k \< 0.7\\, the PSIS estimate and the
  corresponding Monte Carlo standard error estimate are not reliable,
  but increasing the (effective) sample size \\S\\ above 2200 may help
  (this will increase the sample size specific threshold \\(1 - 1 /
  \log\_{10}(2200) \> 0.7\\ and then the bias specific threshold 0.7
  dominates).

- If \\0.7 \<= k \< 1\\, the PSIS estimate and the corresponding Monte
  Carlo standard error have large bias and are not reliable. Increasing
  the sample size may reduce the variability in the \\k\\ estimate,
  which may also result in a lower \\k\\ estimate.

- If \\k \geq 1\\, the target distribution is estimated to have
  non-finite mean. The PSIS estimate and the corresponding Monte Carlo
  standard error are not well defined. Increasing the sample size may
  reduce the variability in \\k\\ estimate, which may also result in a
  lower \\k\\ estimate.

Pareto \\k\\ is also useful as a measure of influence of an observation.
Highly influential observations have high \\k\\ values. Very high \\k\\
values often indicate model misspecification, outliers or mistakes in
data processing. See Section 6 of Gabry et al. (2019) for an example.

### Interpreting `p_loo` when Pareto `k` is large

If \\k \> 0.7\\ then we can also look at the `p_loo` estimate for some
additional information about the problem:

- If `p_loo << p` (the total number of parameters in the model), then
  the model is likely to be misspecified. Posterior predictive checks
  (PPCs) are then likely to also detect the problem. Try using an
  overdispersed model, or add more structural information (nonlinearity,
  mixture model, etc.).

- If `p_loo < p` and the number of parameters `p` is relatively large
  compared to the number of observations (e.g., `p>N/5`), it is likely
  that the model is so flexible or the population prior so weak that
  it’s difficult to predict the left out observation (even for the true
  model). This happens, for example, in the simulated 8 schools (in
  VGG2017), random effect models with a few observations per random
  effect, and Gaussian processes and spatial models with short
  correlation lengths.

- If `p_loo > p`, then the model is likely to be badly misspecified. If
  the number of parameters `p<<N`, then PPCs are also likely to detect
  the problem. See the case study at
  <https://avehtari.github.io/modelselection/roaches.html> for an
  example. If `p` is relatively large compared to the number of
  observations, say `p>N/5` (more accurately we should count number of
  observations influencing each parameter as in hierarchical models some
  groups may have few observations and other groups many), it is
  possible that PPCs won't detect the problem.

## elpd_diff

`elpd_diff` is the difference in `elpd_loo` for two models. If more than
two models are compared, the difference is computed relative to the
model with highest `elpd_loo`.

## se_diff

The standard error of component-wise differences of elpd_loo (Eq 24 in
VGG2017) between two models. This SE is *smaller* than the SE for
individual models due to correlation (i.e., if some observations are
easier and some more difficult to predict for all models).

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
[doi:10.1214/25-BA1569](https://doi.org/10.1214/25-BA1569) .

Gabry, J. , Simpson, D. , Vehtari, A. , Betancourt, M. and Gelman, A.
(2019), Visualization in Bayesian workflow. *J. R. Stat. Soc. A*, 182:
389-402. doi:10.1111/rssa.12378 ([journal
version](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssa.12378),
[preprint arXiv:1709.01449](https://arxiv.org/abs/1709.01449), [code on
GitHub](https://github.com/jgabry/bayes-vis-paper))
