# loo 0.1.4
* Don't call functions from **parallel** package if `cores=1` (This should make
it easier for other package authors using **loo** to write tests that run the
loo function.)
* Fix bug causing `psislw` to return summary statistic rather than smoothed
weights vector. This bug only occurred if `psislw` was called directly and was
not a problem when `psislw` was used internally to compute LOO.

# loo 0.1.3
This update provides several important improvements, most notably an alternative
method for specifying the pointwise log-likelihood that reduces memory usage 
and allows for **loo** to be used with larger datasets. This update also makes
it easier to to incorporate **loo**'s functionality into other packages.

* Add Ben Goodrich as contributor
* S3 generics and `matrix` and `function` methods for both `loo` and `waic`. 
The matrix method provide the same functionality as in previous versions of 
**loo** (taking a log-likelihood matrix as the input). The function method 
allows the user to provide a function for computing the log-likelihood from 
the data and posterior draws (which are also provided by the user). The function
method is less memory intensive and should make it possible to use **loo** for 
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
