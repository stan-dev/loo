# loo 2.6.0.9000

### New features

* `E_loo` now allows `type="sd"`. 


### Bug fixes

* Fix bug in `E_loo` when `type=variance`. 
 
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
