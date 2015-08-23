# loo 0.1.3
This update makes it easier to incorporate **loo**'s functionality into other packages
and also incorporates several changes that reduce memory usage.

* Add Ben Goodrich as contributor
* S3 generics and `matrix` and `function` methods for both `loo` and `waic`. 
The matrix methods provide the same functionality as previous versions of 
**loo** and take a log-likelihood matrix as the input. The function method 
allows the user to provide a function for computing the log-likelihood from 
the data and posterior draws (which are also provided by the user). The function
method is less memory intensive and should make it possible to use **loo** for 
models fit to larger amounts of data than before.
* Separate `plot` and `print` methods. `plot` also provides `label_points` argument,
which, if `TRUE`, will label any Pareto k points greater than 1/2 by the index number
of the corresponding observation.

# loo 0.1.2 
* Updates names in package to reflect name changes in the accompanying 
paper.

# loo 0.1.1
* Better handling of special cases
* Deprecates `loo_and_waic` function in favor of separate functions `loo` and
`waic`
* Deprecates `loo_and_waic_diff`. Use `compare` instead. 

# loo 0.1.0
* Initial release
