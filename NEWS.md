# loo (development version)

* `elpd()`, `crps()`, `scrps()`, `loo_crps()`, `loo_scrps()`, and
  `loo_predictive_metric()` are deprecated in favour of the `measure_*()`
  functions and the `*_pred_measure()` workflow. See `vignette("migration-guide")`
  for a full mapping table.
* New predictive performance API: `insample_pred_measure()`, `loo_pred_measure()`,
  `kfold_pred_measure()`, `test_pred_measure()`, and `pred_measure()` with
  built-in measures via `measure_*()` and [supported_measures_list()].
* `loo_compare()` is now called `model_compare()`. `loo_compare()` remains a
  working alias and is still a generic, so existing code and methods registered
  by other packages keep working.
* `model_compare()` supports every `pred_measure` result --- from
  `loo_pred_measure()`, `kfold_pred_measure()`, `test_pred_measure()`, and
  `insample_pred_measure()` --- with paired differences for all measures common
  to the compared models, optional `rank_by` ranking, utility-scale sign
  conversion for loss measures, and `print(compare, measures = ...)` for
  multi-measure tables by @florence-bockting in #380.
  Measures are matched on bare names, so the source suffix (`_loo`, `_kfold`,
  `_test`, or none for in-sample) is handled transparently. All models in one
  call must share an evaluation source; mixing them is an error. `diag_elpd` is
  reported only for `loo_pred_measure()` comparisons, k-fold comparisons warn on
  differing `K`, and in-sample comparisons warn that they are optimistically
  biased.
* In `model_compare()`, when `rank_by` is not supplied each measure is now
  compared against *its own* best model, so e.g. `mse_diff` can be relative to a
  different model than `elpd_diff`. Rows are still ordered by `"elpd"`, and each
  `{measure}_diff` column has exactly one `0` entry, at that measure's best
  model. Supplying `rank_by` keeps the previous behaviour of pinning the
  top-ranked model as the single reference for every measure. The reference used
  per measure is recorded in the `compare_reference` attribute and shown by
  `print(x, measures = "all")`. Each printed measure table is now sorted by its
  own difference, so the best model on that measure is always the first row and
  the differences run in decreasing order.
* `print()` on a `pred_measure` comparison now always names the reference model
  each set of differences was computed against, including the default case where
  every measure keeps its own best model. The line is written to standard output
  rather than emitted as a message, so it survives `suppressMessages()` and
  knitr chunks with `message = FALSE`. For `loo_pred_measure()` comparisons the
  PSIS-LOO Pareto k diagnostic moved out of the per-measure difference tables
  into a block above them: it describes a model's LOO approximation, not any one
  measure or pair of models, so it is now reported once per model. The
  `diag_elpd` column on the returned object is unchanged, and `print()` for
  plain `"loo"` comparisons is unchanged. Printed sentences are wrapped so that
  no line of output exceeds 80 characters.
* Bug fix: the `diag_diff` column tested `elpd_diff > -4` rather than
  `|elpd_diff| < 4`, so every positive difference was flagged as a small one.
  This was only reachable with `rank_by` set to a model that is not the best
  model, where differences against the reference can be positive.
* In `model_compare()`, `rank_by` now also accepts a **model name** (as shown in
  the `model` column). The named model then becomes the single reference model
  for every measure, whether or not it is the best one, while rows stay ordered
  by `"elpd"`. This form also works for plain `"loo"` comparisons, where
  `elpd_diff` is then relative to the named model. The pinned model is recorded
  in the `rank_by` attribute; a name matching both a measure and a model is
  treated as the measure, with a warning.
* The `diff_method` value `"pairwise"` is now called `"measure_specific"`. Every
  model comparison is pairwise, and the `"sum"`/`"mean"` methods are the ones
  built from paired pointwise differences, so the old name described the one
  category that is *not* a paired pointwise standard deviation. `diff_method` is
  read-only metadata in `measure_info`, so no user code that passes
  arguments is affected; `custom_se_fn` still accepts only `"sum"` and `"mean"`.
* Bug fix: in `model_compare()`, `mlpd` was registered with
  `diff_method = "sum"` although its estimate is the *mean* of the pointwise log
  predictive densities. `mlpd_diff` was therefore reported as the sum of the
  pointwise differences (a factor of `N` too large, and inconsistent with the
  difference of the reported `mlpd` estimates), with a matching `mlpd_se_diff`.
  `mlpd` now uses `diff_method = "mean"`.
* `model_compare()` gains a `custom_se_fn` argument controlling how the standard
  error of a difference is computed for a **custom** measure. It accepts a
  function called as `custom_se_fn(ref, cmp)`, the shorthands `"sum"` and
  `"mean"` for the paired pointwise formulas, or `NULL` to report the difference
  with an `NA` standard error; with two or more custom measures, pass a list
  named by measure. It is required whenever a custom measure is compared ---
  nothing is inferred from a measure's values any more, and custom measures now
  carry `diff_method = "custom"` instead of `"auto"`. The previous
  `attr(my_fun, "se_diff_fun")` route has been removed.
* A custom measure can declare that it is a loss with
  `attr(my_fun, "measure_loss") <- TRUE`, alongside `attr(my_fun,
  "measure_name")`. `model_compare()` then flips its differences onto the
  utility scale and orders `rank_by` by the lowest loss, as it does for built-in
  loss measures; without the declaration a custom measure is still treated as a
  utility. The declaration is recorded as `loss` in the `measure_info`
  attribute, and models that disagree on it cannot be compared.
* Measure values are always stored on the measure's own scale; `loss` is the
  single setting describing a measure's orientation. Measures no longer take a
  `higher_is_better` argument, `control` no longer accepts one, and results no
  longer carry a `measure_higher_is_better` attribute. Comparisons were already
  invariant to the setting --- `model_compare()` reports every measure on the
  utility scale regardless --- so no comparison output changes. Negate the
  values yourself if you want a loss stored as a utility.
* Bug fix: `rps` was registered as a utility and `srps` as a loss, but
  `measure_rps()` returns the Gneiting & Raftery (2007) unscaled score, where
  lower is better, and `measure_srps()` the Bolin & Wallin (2023) scaled score,
  where higher is better. `model_compare()` therefore flipped `srps_diff` and
  left `rps_diff` unflipped, and `rank_by = "rps"` ordered models backwards.
  The registry now records `rps` as a loss and `srps` as a utility.
* The per-measure metadata attribute `measure_compare_meta` is now called
  `measure_info`.
* `control` entries now work for custom measures instead of failing with
  `object 'measure_<name>' of mode 'function' was not found`. A `control` entry
  naming no measure being computed now warns instead of aborting.
* `print()` on a multi-measure comparison now labels each measure's standard
  error column `{measure}_se_diff`, matching the column name on the returned
  data frame, instead of relabelling it `se_diff`. ELPD tables are unchanged.
* Update user messages in `print()` by @ishaan-arora-1, @florence-bockting in 
#328.

# loo 2.10.1

* Revert behavior of internal `psis_smooth_tail()` function that led to several
failing tests not caught by CRAN reverse dependency checks. This change should
not directly affect individual users. by @VisruthSK in #381
* `print.compare.loo()` regains a `simplify` argument for showing the full
comparison table with `simplify = FALSE` by @florence-bockting in #383.

# loo 2.10.0

* Updates to `loo_compare` output by @jgabry, @avehtari, @florence-bockting in #300: 
  - new output columns `p_worse`, `diag_diff`, `diag_elpd`
  - returned object is now a `data.frame` instead of a `matrix`
* Add new `kfold.print` method by @florence-bockting in #342
* Add touchstone performance benchmarks by @VisruthSK in #352
* Use `posterior::gpdfit` and `posterior::qgeneralized_pareto` by @avehtari in #305
* Fix bug in IS method error message + documentation typos by @ishaan-arora-1 in #326
* Shortened `?loo_subsample` title and improved description by @vinniott in #339
* Add AI contribution policy by @florence-bockting in #361


# loo 2.9.0

* Avoid under and overflows in stacking by @avehtari in #273
* Fix `kfold_split_stratified()` when a group has 1 observation by @jgabry in #278
* Fix `plot_diagnostic(..., label_points = TRUE)` by @annariha in #288
* Use testthat 3e and rely on posterior for ESS by @VisruthSK in #289
* Fixed NAs in `psis_r_eff()`  by @VisruthSK in #301
* Print unequal sample sizes in loo_compare by @VisruthSK in #307
* New website theme and pkgdown workflow by @VisruthSK in #292
* Small fix in `loo_model_weights()` doc by @avehtari in #276
* loo_moment_match.R: fix doc for default k_threshold by @jgabry in #279
* Update scrps reference and improve doc by @avehtari in #280
* Fix url in vignette by @jgabry in #282
* Update stacking citation in inst/CITATION by @jgabry in #284
* Added contribution section. by @VisruthSK in #286
* Update LOO uncertainty paper to use BA doi by @avehtari in #311
* Update documentation for `E_loo()` function by @avehtari in #312


# loo 2.8.0

* make E_loo Pareto-k diagnostic more robust by @avehtari in #251
* update psis paper reference by @avehtari in #252
* update PSIS references in vignettes by @jgabry in #254
* fix loo_moment_match p_loo computation by @avehtari in #257
* fix loo_moment_matching NaN issue by @avehtari in #259
* catch Stan log_prob exceptions inside moment matching by @avehtari in #262
* Fix E_loo_khat error when posterior::pareto_khat returns NA by @jgabry in #264
* update psis ref + some minor typo fixes by @avehtari in #266
* update PSIS ref + link to Nabiximols study for Jacobian correction by @avehtari in #267
* Fix issue with pareto_khat output no longer being a list by @n-kall in #269
* fix equations in loo-glossary by @avehtari in #268

# loo 2.7.0

### Major changes

* __New sample size specific diagnostic threshold for Pareto `k`__. The pre-2022 version
of the [PSIS paper](https://arxiv.org/abs/1507.02646) recommended diagnostic
thresholds of 
`k < 0.5 "good"`, `0.5 <= k < 0.7 "ok"`, 
`0.7 <= k < 1 "bad"`, `k>=1 "very bad"`. 
The 2022 revision of the PSIS paper now recommends 
`k < min(1 - 1/log10(S), 0.7) "good"`, `min(1 - 1/log10(S), 0.7) <= k < 1 "bad"`, 
`k > 1 "very bad"`, where `S` is the sample size. 
There is now one fewer diagnostic threshold (`"ok"` has been removed), and the
most important threshold now depends on the sample size `S`. With sample sizes
`100`, `320`, `1000`, `2200`, `10000` the sample size specific part 
`1 - 1/log10(S)` corresponds to thresholds of `0.5`, `0.6`, `0.67`, `0.7`, `0.75`.
Even if the sample size grows, the bias in the PSIS estimate dominates if 
`0.7 <= k < 1`, and thus the diagnostic threshold for good is capped at 
`0.7` (if `k > 1`, the mean does not exist and bias is not a valid measure). 
The new recommended thresholds are based on more careful bias-variance analysis
of PSIS based on truncated Pareto sums theory. For those who use the Stan
default 4000 posterior draws, the `0.7` threshold will be roughly the same, but
there will be fewer warnings as there will be no diagnostic message for `0.5 <=
k < 0.7`. Those who use smaller sample sizes may see diagnostic messages with a
threshold less than `0.7`, and they can simply increase the sample size to about
`2200` to get the threshold to `0.7`. 

* __No more warnings if the `r_eff` argument is not provided__, and the
default is now `r_eff = 1`. The summary print output showing MCSE and ESS now
shows diagnostic information on the range of `r_eff`. The change was made to
reduce unnecessary warnings. The use of `r_eff` does not change the expected
value of `elpd_loo`, `p_loo`, and Pareto `k`, and is needed only to estimate
MCSE and ESS. Thus it is better to show the diagnostic information about `r_eff`
only when MCSE and ESS values are shown.

### Other changes

* Make Pareto `k` Inf if it is NA by @topipa in #224
* Fix bug in `E_loo()` when type is variance by @jgabry in #226
* `E_loo()` now allows `type="sd"` by @jgabry in #226
* update array syntax in vignettes by @jgabry in #229
* Fix unbalanced knitr backticks by @jgabry in #232
* include cc-by 4.0 license for documentation  by @jgabry in #216
* Add order statistic warning by @yannmclatchie in #230
* `pointwise()` convenience function for extracting pointwise estimates by @jgabry in #241
* use new `k` threshold by @avehtari in #235
* simplify `mcse_elpd` using log-normal approximation by @avehtari in #246
* show NA for `n_eff/ESS` if `k > k_threshold` by @avehtari in #248
* improved `E_loo()` Pareto-k diagnostics by @avehtari in #247
* Doc improvement in `loo_subsample.R` by @avehtari in #238
* Fix typo and deprecations in LFO vignette by @jgabry in #244
* Register internal S3 methods by @jgabry in #239
* Avoid R cmd check NOTEs about some internal functions by @jgabry in #240
* fix R cmd check note due to importance_sampling roxygen template by @jgabry in #233
* fix R cmd check notes by @jgabry in #242


 
# loo 2.6.0

### New features 

* New `loo_predictive_metric()` function for computing estimates of leave-one-out
predictive metrics: mean absolute error, mean squared error and root mean
squared error for continuous predictions, and accuracy and balanced accuracy for
binary classification. (#202, @LeeviLindgren)

* New functions `crps()`, `scrps()`, `loo_crps()`, and `loo_scrps()` for
computing the (scaled) continuously ranked probability score. (#203, @LeeviLindgren)

* New vignette "Mixture IS leave-one-out cross-validation for high-dimensional Bayesian models." This is a demonstration of the mixture estimators proposed by [Silva and Zanella (2022)](https://arxiv.org/abs/2209.09190). (#210)

### Bug fixes

* Minor fix to model names displayed by `loo_model_weights()` to make them consistent with `loo_compare()`. (#217)


# loo 2.5.1

* Fix R CMD check error on M1 Mac

# loo 2.5.0

### Improvements

* New [Frequently Asked Questions page](https://mc-stan.org/loo/articles/online-only/faq.html) on the package website. (#143)

* Speed improvement from simplifying the normalization when fitting the 
generalized Pareto distribution. (#187, @sethaxen)

* Added parallel likelihood computation to speedup `loo_subsample()` when using posterior approximations. (#171, @kdubovikov)

* Switch unit tests from Travis to GitHub Actions. (#164)

### Bug fixes 

* Fixed a bug causing the normalizing constant of the PSIS (log) weights not 
to get updated when performing moment matching with `save_psis = TRUE` (#166, @fweber144).

# loo 2.4.1

* Fixed issue reported by CRAN where one of the vignettes errored on an M1 Mac
due to RStan's dependency on V8. 

# loo 2.4.0

### Bug fixes 

* Fixed a bug in `relative_eff.function()` that caused an error on Windows when
using multiple cores. (#152) 

* Fixed a potential numerical issue in `loo_moment_match()` with `split=TRUE`. (#153)

* Fixed potential integer overflow with `loo_moment_match()`. (#155, @ecmerkle)

* Fixed `relative_eff()` when used with a `posterior::draws_array`. (#161, @rok-cesnovar)

### New features 

* New generic function `elpd()` (and methods for matrices and arrays) for
computing expected log predictive density of new data or log predictive density
of observed data. A new vignette demonstrates using this function when doing
K-fold CV with rstan. (#159, @bnicenboim)


# loo 2.3.1

* Fixed a bug in `loo_moment_match()` that prevented `...` arguments from being
used correctly. (#149)


# loo 2.3.0

* Added Topi Paananen and Paul Bürkner as coauthors.

* New function `loo_moment_match()` (and new vignette), which can be used to
update a `loo` object when Pareto k estimates are large. (#130)

* The log weights provided by the importance sampling functions `psis()`,
`tis()`, and `sis()` no longer have the largest log ratio subtracted from them
when returned to the user. This should be less confusing for anyone using 
the `weights()` method to make an importance sampler. (#112, #146)

* MCSE calculation is now deterministic (#116, #147)


# loo 2.2.0

* Added Mans Magnusson as a coauthor.

* New functions `loo_subsample()` and `loo_approximate_posterior()` (and new
vignette) for doing PSIS-LOO with large data. (#113)

* Added support for standard importance sampling and truncated importance
sampling (functions `sis()` and `tis()`). (#125)

* `compare()` now throws a deprecation warning suggesting `loo_compare()`. (#93)

* A smaller threshold is used when checking the uniqueness of tail values. (#124)

* For WAIC, warnings are only thrown when running `waic()` and not when printing
a `waic` object. (#117, @mcol)

* Use markdown syntax in roxygen documentation wherever possible. (#108)


# loo 2.1.0

* New function `loo_compare()` for model comparison that will eventually replace
the existing `compare()` function. (#93)

* New vignette on LOO for non-factorizable joint Gaussian models. (#75)

* New vignette on "leave-future-out" cross-validation for time series models. (#90)

* New glossary page (use `help("loo-glossary")`) with definitions of key terms. (#81)

* New `se_diff` column in model comparison results. (#78)

* Improved stability of `psis()` when `log_ratios` are very small. (#74)

* Allow `r_eff=NA` to suppress warning when specifying `r_eff` is not applicable 
(i.e., draws not from MCMC). (#72)

* Update effective sample size calculations to match RStan's version. (#85)

* Naming of k-fold helper functions now matches scikit-learn. (#96)

# loo 2.0.0

This is a major release with many changes. Whenever possible we have opted to
deprecate rather than remove old functionality, but it is possible that old code
that accesses elements inside loo objects by position rather than name may
error.

* New package documentation website http://mc-stan.org/loo/ with vignettes, 
function reference, news.

* Updated existing vignette and added two new vignettes demonstrating how to use
the package.

* New function `psis()` replaces `psislw()` (now deprecated). This version 
implements the improvements to the PSIS algorithm described in the latest
version of https://arxiv.org/abs/1507.02646. Additional diagnostic
information is now also provided, including PSIS effective sample sizes.

* New `weights()` method for extracting smoothed weights from a `psis` object.
Arguments `log` and `normalize` control whether the weights are returned on the
log scale and whether they are normalized.

* Updated the interface for the `loo()` methods to integrate nicely with the new
PSIS algorithm. Methods for log-likelihood arrays, matrices, and functions 
are provided. Several arguments have changed, particularly for the
`loo.function` method. The documentation at `help("loo")` has been updated to
describe the new behavior.

* The structure of the objects returned by the `loo()` function has also changed
slightly, as described in the __Value__ section at `help("loo", package = "loo")`.

* New function `loo_model_weights()` computes weights for model averaging as
described in https://arxiv.org/abs/1704.02030. Implemented methods include
stacking of predictive distributions, pseudo-BMA weighting or pseudo-BMA+
weighting with the Bayesian bootstrap.

* Setting `options(loo.cores=...)` is now deprecated in favor of
`options(mc.cores=...)`. For now, if both the `loo.cores` and `mc.cores` options
have been set, preference will be given to `loo.cores` until it is removed in a
future release. (thanks to @cfhammill)

* New functions `example_loglik_array()` and `example_loglik_matrix()` that
provide objects to use in examples and tests.

* When comparing more than two models with `compare()`, the first column of the
output is now the `elpd` difference from the model in the first row.

* New helper functions for splitting observations for K-fold CV:
`kfold_split_random()`, `kfold_split_balanced()`, `kfold_split_stratified()`.
Additional helper functions for implementing K-fold CV will be included in
future releases.


# loo 1.1.0

* Introduce the `E_loo()` function for computing weighted expectations (means,
variances, quantiles).

# loo 1.0.0

* `pareto_k_table()` and `pareto_k_ids()` convenience functions for quickly
identifying problematic observations
* pareto k values now grouped into `(-Inf, 0.5]`, `(0.5, 0.7]`, `(0.7, 1]`, 
`(1, Inf)` (didn't used to include 0.7)
* warning messages are now issued by `psislw()` instead of `print.loo`
* `print.loo()` shows a table of pareto k estimates (if any k > 0.7)
* Add argument to `compare()` to allow loo objects to be provided in a list
rather than in `'...'`
* Update references to point to published paper

# loo 0.1.6

* GitHub repository moved from @jgabry to @stan-dev
* Better error messages from `extract_log_lik()`
* Fix example code in vignette (thanks to GitHub user @krz)

# loo 0.1.5

* Add warnings if any p_waic estimates are greather than 0.4
* Improve line coverage of tests to 100%
* Update references in documentation
* Remove model weights from `compare()`. In previous versions of __loo__ model
weights were also reported by `compare()`. We have removed the weights because
they were based only on the point estimate of the elpd values ignoring the
uncertainty. We are currently working on something similar to these weights that
also accounts for uncertainty, which will be included in future versions of
__loo__.

# loo 0.1.4

This update makes it easier for other package authors using __loo__ to write
tests that involve running the `loo` function. It also includes minor bug
fixes and additional unit tests. Highlights:

* Don't call functions from __parallel__ package if `cores=1`.
* Return entire vector/matrix of smoothed weights rather than a summary
statistic when `psislw` function is called in an interactive session.
* Test coverage > 80%

# loo 0.1.3

This update provides several important improvements, most notably an alternative
method for specifying the pointwise log-likelihood that reduces memory usage 
and allows for __loo__ to be used with larger datasets. This update also makes
it easier to to incorporate __loo__'s functionality into other packages.

* Add Ben Goodrich as contributor
* S3 generics and `matrix` and `function` methods for both `loo()` and `waic()`. 
The matrix method provide the same functionality as in previous versions of 
__loo__ (taking a log-likelihood matrix as the input). The function method 
allows the user to provide a function for computing the log-likelihood from 
the data and posterior draws (which are also provided by the user). The function
method is less memory intensive and should make it possible to use __loo__ for 
models fit to larger amounts of data than before.
* Separate `plot` and `print` methods. `plot` also provides `label_points` 
argument, which, if `TRUE`, will label any Pareto `k` points greater than 
1/2 by the index number of the corresponding observation. The plot method 
also now warns about `Inf`/`NA`/`NaN` values of `k` that are not shown in 
the plot. 
* `compare` now returns model weights and accepts more than two inputs.
* Allow setting number of cores using `options(loo.cores = NUMBER)`. 

# loo 0.1.2 

* Updates names in package to reflect name changes in the accompanying paper.

# loo 0.1.1

* Better handling of special cases
* Deprecates `loo_and_waic` function in favor of separate functions `loo` and
`waic`
* Deprecates `loo_and_waic_diff`. Use `compare` instead. 

# loo 0.1.0

* Initial release
