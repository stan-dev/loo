# Changelog

## loo 2.8.0.9000

Items for next release go here

## loo 2.8.0

CRAN release: 2024-07-03

- make E_loo Pareto-k diagnostic more robust by
  [@avehtari](https://github.com/avehtari) in
  [\#251](https://github.com/stan-dev/loo/issues/251)
- update psis paper reference by
  [@avehtari](https://github.com/avehtari) in
  [\#252](https://github.com/stan-dev/loo/issues/252)
- update PSIS references in vignettes by
  [@jgabry](https://github.com/jgabry) in
  [\#254](https://github.com/stan-dev/loo/issues/254)
- fix loo_moment_match p_loo computation by
  [@avehtari](https://github.com/avehtari) in
  [\#257](https://github.com/stan-dev/loo/issues/257)
- fix loo_moment_matching NaN issue by
  [@avehtari](https://github.com/avehtari) in
  [\#259](https://github.com/stan-dev/loo/issues/259)
- catch Stan log_prob exceptions inside moment matching by
  [@avehtari](https://github.com/avehtari) in
  [\#262](https://github.com/stan-dev/loo/issues/262)
- Fix E_loo_khat error when posterior::pareto_khat returns NA by
  [@jgabry](https://github.com/jgabry) in
  [\#264](https://github.com/stan-dev/loo/issues/264)
- update psis ref + some minor typo fixes by
  [@avehtari](https://github.com/avehtari) in
  [\#266](https://github.com/stan-dev/loo/issues/266)
- update PSIS ref + link to Nabiximols study for Jacobian correction by
  [@avehtari](https://github.com/avehtari) in
  [\#267](https://github.com/stan-dev/loo/issues/267)
- Fix issue with pareto_khat output no longer being a list by
  [@n-kall](https://github.com/n-kall) in
  [\#269](https://github.com/stan-dev/loo/issues/269)
- fix equations in loo-glossary by
  [@avehtari](https://github.com/avehtari) in
  [\#268](https://github.com/stan-dev/loo/issues/268)

## loo 2.7.0

CRAN release: 2024-02-24

#### Major changes

- **New sample size specific diagnostic threshold for Pareto `k`**. The
  pre-2022 version of the [PSIS paper](https://arxiv.org/abs/1507.02646)
  recommended diagnostic thresholds of `k < 0.5 "good"`,
  `0.5 <= k < 0.7 "ok"`, `0.7 <= k < 1 "bad"`, `k>=1 "very bad"`. The
  2022 revision of the PSIS paper now recommends
  `k < min(1 - 1/log10(S), 0.7) "good"`,
  `min(1 - 1/log10(S), 0.7) <= k < 1 "bad"`, `k > 1 "very bad"`, where
  `S` is the sample size. There is now one fewer diagnostic threshold
  (`"ok"` has been removed), and the most important threshold now
  depends on the sample size `S`. With sample sizes `100`, `320`,
  `1000`, `2200`, `10000` the sample size specific part `1 - 1/log10(S)`
  corresponds to thresholds of `0.5`, `0.6`, `0.67`, `0.7`, `0.75`. Even
  if the sample size grows, the bias in the PSIS estimate dominates if
  `0.7 <= k < 1`, and thus the diagnostic threshold for good is capped
  at `0.7` (if `k > 1`, the mean does not exist and bias is not a valid
  measure). The new recommended thresholds are based on more careful
  bias-variance analysis of PSIS based on truncated Pareto sums theory.
  For those who use the Stan default 4000 posterior draws, the `0.7`
  threshold will be roughly the same, but there will be fewer warnings
  as there will be no diagnostic message for `0.5 <= k < 0.7`. Those who
  use smaller sample sizes may see diagnostic messages with a threshold
  less than `0.7`, and they can simply increase the sample size to about
  `2200` to get the threshold to `0.7`.

- **No more warnings if the `r_eff` argument is not provided**, and the
  default is now `r_eff = 1`. The summary print output showing MCSE and
  ESS now shows diagnostic information on the range of `r_eff`. The
  change was made to reduce unnecessary warnings. The use of `r_eff`
  does not change the expected value of `elpd_loo`, `p_loo`, and Pareto
  `k`, and is needed only to estimate MCSE and ESS. Thus it is better to
  show the diagnostic information about `r_eff` only when MCSE and ESS
  values are shown.

#### Other changes

- Make Pareto `k` Inf if it is NA by
  [@topipa](https://github.com/topipa) in
  [\#224](https://github.com/stan-dev/loo/issues/224)
- Fix bug in [`E_loo()`](https://mc-stan.org/loo/dev/reference/E_loo.md)
  when type is variance by [@jgabry](https://github.com/jgabry) in
  [\#226](https://github.com/stan-dev/loo/issues/226)
- [`E_loo()`](https://mc-stan.org/loo/dev/reference/E_loo.md) now allows
  `type="sd"` by [@jgabry](https://github.com/jgabry) in
  [\#226](https://github.com/stan-dev/loo/issues/226)
- update array syntax in vignettes by
  [@jgabry](https://github.com/jgabry) in
  [\#229](https://github.com/stan-dev/loo/issues/229)
- Fix unbalanced knitr backticks by [@jgabry](https://github.com/jgabry)
  in [\#232](https://github.com/stan-dev/loo/issues/232)
- include cc-by 4.0 license for documentation by
  [@jgabry](https://github.com/jgabry) in
  [\#216](https://github.com/stan-dev/loo/issues/216)
- Add order statistic warning by
  [@yannmclatchie](https://github.com/yannmclatchie) in
  [\#230](https://github.com/stan-dev/loo/issues/230)
- [`pointwise()`](https://mc-stan.org/loo/dev/reference/pointwise.md)
  convenience function for extracting pointwise estimates by
  [@jgabry](https://github.com/jgabry) in
  [\#241](https://github.com/stan-dev/loo/issues/241)
- use new `k` threshold by [@avehtari](https://github.com/avehtari) in
  [\#235](https://github.com/stan-dev/loo/issues/235)
- simplify `mcse_elpd` using log-normal approximation by
  [@avehtari](https://github.com/avehtari) in
  [\#246](https://github.com/stan-dev/loo/issues/246)
- show NA for `n_eff/ESS` if `k > k_threshold` by
  [@avehtari](https://github.com/avehtari) in
  [\#248](https://github.com/stan-dev/loo/issues/248)
- improved [`E_loo()`](https://mc-stan.org/loo/dev/reference/E_loo.md)
  Pareto-k diagnostics by [@avehtari](https://github.com/avehtari) in
  [\#247](https://github.com/stan-dev/loo/issues/247)
- Doc improvement in `loo_subsample.R` by
  [@avehtari](https://github.com/avehtari) in
  [\#238](https://github.com/stan-dev/loo/issues/238)
- Fix typo and deprecations in LFO vignette by
  [@jgabry](https://github.com/jgabry) in
  [\#244](https://github.com/stan-dev/loo/issues/244)
- Register internal S3 methods by [@jgabry](https://github.com/jgabry)
  in [\#239](https://github.com/stan-dev/loo/issues/239)
- Avoid R cmd check NOTEs about some internal functions by
  [@jgabry](https://github.com/jgabry) in
  [\#240](https://github.com/stan-dev/loo/issues/240)
- fix R cmd check note due to importance_sampling roxygen template by
  [@jgabry](https://github.com/jgabry) in
  [\#233](https://github.com/stan-dev/loo/issues/233)
- fix R cmd check notes by [@jgabry](https://github.com/jgabry) in
  [\#242](https://github.com/stan-dev/loo/issues/242)

## loo 2.6.0

CRAN release: 2023-03-31

#### New features

- New
  [`loo_predictive_metric()`](https://mc-stan.org/loo/dev/reference/loo_predictive_metric.md)
  function for computing estimates of leave-one-out predictive metrics:
  mean absolute error, mean squared error and root mean squared error
  for continuous predictions, and accuracy and balanced accuracy for
  binary classification.
  ([\#202](https://github.com/stan-dev/loo/issues/202),
  [@LeeviLindgren](https://github.com/LeeviLindgren))

- New functions
  [`crps()`](https://mc-stan.org/loo/dev/reference/crps.md),
  [`scrps()`](https://mc-stan.org/loo/dev/reference/crps.md),
  [`loo_crps()`](https://mc-stan.org/loo/dev/reference/crps.md), and
  [`loo_scrps()`](https://mc-stan.org/loo/dev/reference/crps.md) for
  computing the (scaled) continuously ranked probability score.
  ([\#203](https://github.com/stan-dev/loo/issues/203),
  [@LeeviLindgren](https://github.com/LeeviLindgren))

- New vignette “Mixture IS leave-one-out cross-validation for
  high-dimensional Bayesian models.” This is a demonstration of the
  mixture estimators proposed by [Silva and Zanella
  (2022)](https://arxiv.org/abs/2209.09190).
  ([\#210](https://github.com/stan-dev/loo/issues/210))

#### Bug fixes

- Minor fix to model names displayed by
  [`loo_model_weights()`](https://mc-stan.org/loo/dev/reference/loo_model_weights.md)
  to make them consistent with
  [`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md).
  ([\#217](https://github.com/stan-dev/loo/issues/217))

## loo 2.5.1

CRAN release: 2022-03-24

- Fix R CMD check error on M1 Mac

## loo 2.5.0

CRAN release: 2022-03-16

#### Improvements

- New [Frequently Asked Questions
  page](https://mc-stan.org/loo/articles/online-only/faq.html) on the
  package website. ([\#143](https://github.com/stan-dev/loo/issues/143))

- Speed improvement from simplifying the normalization when fitting the
  generalized Pareto distribution.
  ([\#187](https://github.com/stan-dev/loo/issues/187),
  [@sethaxen](https://github.com/sethaxen))

- Added parallel likelihood computation to speedup
  [`loo_subsample()`](https://mc-stan.org/loo/dev/reference/loo_subsample.md)
  when using posterior approximations.
  ([\#171](https://github.com/stan-dev/loo/issues/171),
  [@kdubovikov](https://github.com/kdubovikov))

- Switch unit tests from Travis to GitHub Actions.
  ([\#164](https://github.com/stan-dev/loo/issues/164))

#### Bug fixes

- Fixed a bug causing the normalizing constant of the PSIS (log) weights
  not to get updated when performing moment matching with
  `save_psis = TRUE`
  ([\#166](https://github.com/stan-dev/loo/issues/166),
  [@fweber144](https://github.com/fweber144)).

## loo 2.4.1

CRAN release: 2020-12-09

- Fixed issue reported by CRAN where one of the vignettes errored on an
  M1 Mac due to RStan’s dependency on V8.

## loo 2.4.0

CRAN release: 2020-12-05

#### Bug fixes

- Fixed a bug in
  [`relative_eff.function()`](https://mc-stan.org/loo/dev/reference/relative_eff.md)
  that caused an error on Windows when using multiple cores.
  ([\#152](https://github.com/stan-dev/loo/issues/152))

- Fixed a potential numerical issue in
  [`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
  with `split=TRUE`.
  ([\#153](https://github.com/stan-dev/loo/issues/153))

- Fixed potential integer overflow with
  [`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md).
  ([\#155](https://github.com/stan-dev/loo/issues/155),
  [@ecmerkle](https://github.com/ecmerkle))

- Fixed
  [`relative_eff()`](https://mc-stan.org/loo/dev/reference/relative_eff.md)
  when used with a
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html).
  ([\#161](https://github.com/stan-dev/loo/issues/161),
  [@rok-cesnovar](https://github.com/rok-cesnovar))

#### New features

- New generic function
  [`elpd()`](https://mc-stan.org/loo/dev/reference/elpd.md) (and methods
  for matrices and arrays) for computing expected log predictive density
  of new data or log predictive density of observed data. A new vignette
  demonstrates using this function when doing K-fold CV with rstan.
  ([\#159](https://github.com/stan-dev/loo/issues/159),
  [@bnicenboim](https://github.com/bnicenboim))

## loo 2.3.1

CRAN release: 2020-07-14

- Fixed a bug in
  [`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
  that prevented `...` arguments from being used correctly.
  ([\#149](https://github.com/stan-dev/loo/issues/149))

## loo 2.3.0

CRAN release: 2020-07-07

- Added Topi Paananen and Paul Bürkner as coauthors.

- New function
  [`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
  (and new vignette), which can be used to update a `loo` object when
  Pareto k estimates are large.
  ([\#130](https://github.com/stan-dev/loo/issues/130))

- The log weights provided by the importance sampling functions
  [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md),
  [`tis()`](https://mc-stan.org/loo/dev/reference/tis.md), and
  [`sis()`](https://mc-stan.org/loo/dev/reference/sis.md) no longer have
  the largest log ratio subtracted from them when returned to the user.
  This should be less confusing for anyone using the
  [`weights()`](https://rdrr.io/r/stats/weights.html) method to make an
  importance sampler.
  ([\#112](https://github.com/stan-dev/loo/issues/112),
  [\#146](https://github.com/stan-dev/loo/issues/146))

- MCSE calculation is now deterministic
  ([\#116](https://github.com/stan-dev/loo/issues/116),
  [\#147](https://github.com/stan-dev/loo/issues/147))

## loo 2.2.0

CRAN release: 2019-12-19

- Added Mans Magnusson as a coauthor.

- New functions
  [`loo_subsample()`](https://mc-stan.org/loo/dev/reference/loo_subsample.md)
  and
  [`loo_approximate_posterior()`](https://mc-stan.org/loo/dev/reference/loo_approximate_posterior.md)
  (and new vignette) for doing PSIS-LOO with large data.
  ([\#113](https://github.com/stan-dev/loo/issues/113))

- Added support for standard importance sampling and truncated
  importance sampling (functions
  [`sis()`](https://mc-stan.org/loo/dev/reference/sis.md) and
  [`tis()`](https://mc-stan.org/loo/dev/reference/tis.md)).
  ([\#125](https://github.com/stan-dev/loo/issues/125))

- [`compare()`](https://mc-stan.org/loo/dev/reference/compare.md) now
  throws a deprecation warning suggesting
  [`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md).
  ([\#93](https://github.com/stan-dev/loo/issues/93))

- A smaller threshold is used when checking the uniqueness of tail
  values. ([\#124](https://github.com/stan-dev/loo/issues/124))

- For WAIC, warnings are only thrown when running
  [`waic()`](https://mc-stan.org/loo/dev/reference/waic.md) and not when
  printing a `waic` object.
  ([\#117](https://github.com/stan-dev/loo/issues/117),
  [@mcol](https://github.com/mcol))

- Use markdown syntax in roxygen documentation wherever possible.
  ([\#108](https://github.com/stan-dev/loo/issues/108))

## loo 2.1.0

CRAN release: 2019-03-13

- New function
  [`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md)
  for model comparison that will eventually replace the existing
  [`compare()`](https://mc-stan.org/loo/dev/reference/compare.md)
  function. ([\#93](https://github.com/stan-dev/loo/issues/93))

- New vignette on LOO for non-factorizable joint Gaussian models.
  ([\#75](https://github.com/stan-dev/loo/issues/75))

- New vignette on “leave-future-out” cross-validation for time series
  models. ([\#90](https://github.com/stan-dev/loo/issues/90))

- New glossary page (use
  [`help("loo-glossary")`](https://mc-stan.org/loo/dev/reference/loo-glossary.md))
  with definitions of key terms.
  ([\#81](https://github.com/stan-dev/loo/issues/81))

- New `se_diff` column in model comparison results.
  ([\#78](https://github.com/stan-dev/loo/issues/78))

- Improved stability of
  [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md) when
  `log_ratios` are very small.
  ([\#74](https://github.com/stan-dev/loo/issues/74))

- Allow `r_eff=NA` to suppress warning when specifying `r_eff` is not
  applicable (i.e., draws not from MCMC).
  ([\#72](https://github.com/stan-dev/loo/issues/72))

- Update effective sample size calculations to match RStan’s version.
  ([\#85](https://github.com/stan-dev/loo/issues/85))

- Naming of k-fold helper functions now matches scikit-learn.
  ([\#96](https://github.com/stan-dev/loo/issues/96))

## loo 2.0.0

CRAN release: 2018-04-11

This is a major release with many changes. Whenever possible we have
opted to deprecate rather than remove old functionality, but it is
possible that old code that accesses elements inside loo objects by
position rather than name may error.

- New package documentation website <http://mc-stan.org/loo/> with
  vignettes, function reference, news.

- Updated existing vignette and added two new vignettes demonstrating
  how to use the package.

- New function [`psis()`](https://mc-stan.org/loo/dev/reference/psis.md)
  replaces [`psislw()`](https://mc-stan.org/loo/dev/reference/psislw.md)
  (now deprecated). This version implements the improvements to the PSIS
  algorithm described in the latest version of
  <https://arxiv.org/abs/1507.02646>. Additional diagnostic information
  is now also provided, including PSIS effective sample sizes.

- New [`weights()`](https://rdrr.io/r/stats/weights.html) method for
  extracting smoothed weights from a `psis` object. Arguments `log` and
  `normalize` control whether the weights are returned on the log scale
  and whether they are normalized.

- Updated the interface for the
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) methods to
  integrate nicely with the new PSIS algorithm. Methods for
  log-likelihood arrays, matrices, and functions are provided. Several
  arguments have changed, particularly for the `loo.function` method.
  The documentation at
  [`help("loo")`](https://mc-stan.org/loo/dev/reference/loo.md) has been
  updated to describe the new behavior.

- The structure of the objects returned by the
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) function has
  also changed slightly, as described in the **Value** section at
  [`help("loo", package = "loo")`](https://mc-stan.org/loo/dev/reference/loo.md).

- New function
  [`loo_model_weights()`](https://mc-stan.org/loo/dev/reference/loo_model_weights.md)
  computes weights for model averaging as described in
  <https://arxiv.org/abs/1704.02030>. Implemented methods include
  stacking of predictive distributions, pseudo-BMA weighting or
  pseudo-BMA+ weighting with the Bayesian bootstrap.

- Setting `options(loo.cores=...)` is now deprecated in favor of
  `options(mc.cores=...)`. For now, if both the `loo.cores` and
  `mc.cores` options have been set, preference will be given to
  `loo.cores` until it is removed in a future release. (thanks to
  [@cfhammill](https://github.com/cfhammill))

- New functions
  [`example_loglik_array()`](https://mc-stan.org/loo/dev/reference/example_loglik_array.md)
  and
  [`example_loglik_matrix()`](https://mc-stan.org/loo/dev/reference/example_loglik_array.md)
  that provide objects to use in examples and tests.

- When comparing more than two models with
  [`compare()`](https://mc-stan.org/loo/dev/reference/compare.md), the
  first column of the output is now the `elpd` difference from the model
  in the first row.

- New helper functions for splitting observations for K-fold CV:
  [`kfold_split_random()`](https://mc-stan.org/loo/dev/reference/kfold-helpers.md),
  `kfold_split_balanced()`,
  [`kfold_split_stratified()`](https://mc-stan.org/loo/dev/reference/kfold-helpers.md).
  Additional helper functions for implementing K-fold CV will be
  included in future releases.

## loo 1.1.0

CRAN release: 2017-03-27

- Introduce the
  [`E_loo()`](https://mc-stan.org/loo/dev/reference/E_loo.md) function
  for computing weighted expectations (means, variances, quantiles).

## loo 1.0.0

CRAN release: 2016-12-16

- [`pareto_k_table()`](https://mc-stan.org/loo/dev/reference/pareto-k-diagnostic.md)
  and
  [`pareto_k_ids()`](https://mc-stan.org/loo/dev/reference/pareto-k-diagnostic.md)
  convenience functions for quickly identifying problematic observations
- pareto k values now grouped into `(-Inf, 0.5]`, `(0.5, 0.7]`,
  `(0.7, 1]`, `(1, Inf)` (didn’t used to include 0.7)
- warning messages are now issued by
  [`psislw()`](https://mc-stan.org/loo/dev/reference/psislw.md) instead
  of `print.loo`
- [`print.loo()`](https://mc-stan.org/loo/dev/reference/print.loo.md)
  shows a table of pareto k estimates (if any k \> 0.7)
- Add argument to
  [`compare()`](https://mc-stan.org/loo/dev/reference/compare.md) to
  allow loo objects to be provided in a list rather than in `'...'`
- Update references to point to published paper

## loo 0.1.6

CRAN release: 2016-03-23

- GitHub repository moved from [@jgabry](https://github.com/jgabry) to
  [@stan-dev](https://github.com/stan-dev)
- Better error messages from
  [`extract_log_lik()`](https://mc-stan.org/loo/dev/reference/extract_log_lik.md)
- Fix example code in vignette (thanks to GitHub user
  [@krz](https://github.com/krz))

## loo 0.1.5

CRAN release: 2016-02-12

- Add warnings if any p_waic estimates are greather than 0.4
- Improve line coverage of tests to 100%
- Update references in documentation
- Remove model weights from
  [`compare()`](https://mc-stan.org/loo/dev/reference/compare.md). In
  previous versions of **loo** model weights were also reported by
  [`compare()`](https://mc-stan.org/loo/dev/reference/compare.md). We
  have removed the weights because they were based only on the point
  estimate of the elpd values ignoring the uncertainty. We are currently
  working on something similar to these weights that also accounts for
  uncertainty, which will be included in future versions of **loo**.

## loo 0.1.4

CRAN release: 2015-12-10

This update makes it easier for other package authors using **loo** to
write tests that involve running the `loo` function. It also includes
minor bug fixes and additional unit tests. Highlights:

- Don’t call functions from **parallel** package if `cores=1`.
- Return entire vector/matrix of smoothed weights rather than a summary
  statistic when `psislw` function is called in an interactive session.
- Test coverage \> 80%

## loo 0.1.3

CRAN release: 2015-09-18

This update provides several important improvements, most notably an
alternative method for specifying the pointwise log-likelihood that
reduces memory usage and allows for **loo** to be used with larger
datasets. This update also makes it easier to to incorporate **loo**’s
functionality into other packages.

- Add Ben Goodrich as contributor
- S3 generics and `matrix` and `function` methods for both
  [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) and
  [`waic()`](https://mc-stan.org/loo/dev/reference/waic.md). The matrix
  method provide the same functionality as in previous versions of
  **loo** (taking a log-likelihood matrix as the input). The function
  method allows the user to provide a function for computing the
  log-likelihood from the data and posterior draws (which are also
  provided by the user). The function method is less memory intensive
  and should make it possible to use **loo** for models fit to larger
  amounts of data than before.
- Separate `plot` and `print` methods. `plot` also provides
  `label_points` argument, which, if `TRUE`, will label any Pareto `k`
  points greater than 1/2 by the index number of the corresponding
  observation. The plot method also now warns about `Inf`/`NA`/`NaN`
  values of `k` that are not shown in the plot.
- `compare` now returns model weights and accepts more than two inputs.
- Allow setting number of cores using `options(loo.cores = NUMBER)`.

## loo 0.1.2

CRAN release: 2015-07-17

- Updates names in package to reflect name changes in the accompanying
  paper.

## loo 0.1.1

- Better handling of special cases
- Deprecates `loo_and_waic` function in favor of separate functions
  `loo` and `waic`
- Deprecates `loo_and_waic_diff`. Use `compare` instead.

## loo 0.1.0

CRAN release: 2015-06-26

- Initial release
