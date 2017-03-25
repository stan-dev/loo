# loo 1.1.0
* Introduce the `E_loo` function for computing weighted expectations (means, variances, quantiles).

# loo 1.0.0
* `pareto_k_table` and `pareto_k_ids` convenience functions for quickly identifying problematic observations
* pareto k values now grouped into `(-Inf, 0.5]`, `(0.5, 0.7]`, `(0.7, 1]`, 
`(1, Inf)` (didn't used to include 0.7)
* warning messages are now issued by `psislw` instead of `print.loo`
* `print.loo` shows a table of pareto k estimates (if any k > 0.7)
* Add argument to `compare` to allow loo objects to be provided in a list rather than in `'...'`
* Update references to point to published paper

# loo 0.1.6
* GitHub repository moved from @jgabry to @stan-dev
* Better error messages from `extract_log_lik`
* Fix example code in vignette (thanks to GitHub user @krz)

# loo 0.1.5
* Add warnings if any p_waic estimates are greather than 0.4
* Improve line coverage of tests to 100%
* Update references in documentation
* Remove model weights from `compare`. <br>
    In previous versions of __loo__ model
weights were also reported by `compare`. We have removed the weights because
they were based only on the point estimate of the elpd values ignoring the
uncertainty. We are currently working on something similar to these weights that
also accounts for uncertainty, which will be included in future versions of
__loo__.

# loo 0.1.4
This update makes it easier for other package authors using __loo__ to write
tests that involve running the `loo` function. It also includes minor bug
fixes and additional unit tests. Highlights:

* Don't call functions from __parallel__ package if `cores=1`.
* Return entire vector/matrix of smoothed weights rather than a summary statistic when `psislw` function is called in an interactive session.
* [Test coverage > 80%](https://codecov.io/github/jgabry/loo?branch=master)

# loo 0.1.3
This update provides several important improvements, most notably an alternative
method for specifying the pointwise log-likelihood that reduces memory usage 
and allows for __loo__ to be used with larger datasets. This update also makes
it easier to to incorporate __loo__'s functionality into other packages.

* Add Ben Goodrich as contributor
* S3 generics and `matrix` and `function` methods for both `loo` and `waic`. 
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
