# loo 0.1.3
This update makes it easier to incorporate **loo**'s functionality into other packages
and also provides a new means of specifying the pointwise log-likelihood that uses
less memory:

* S3 generics and `matrix` and `function` methods for both `loo` and `waic`. 
The matrix methods provide the same functionality as previous versions of 
**loo** and take a log-likelihood matrix as the input. The function method 
allows the user to provide a function for computing the log-likelihood from 
the data and posterior draws (which are also provided by the user). The function
method saves a lot of memory and so should make it possible to use **loo** for 
models fit to large amounts of data.
* Add Ben Goodrich as contributor

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
