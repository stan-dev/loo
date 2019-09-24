#' Generic function for K-fold cross-validation for developers
#'
#' @description For developers of Bayesian modeling packages, **loo** includes
#'   a generic function `kfold()` so that methods may be defined for K-fold
#'   CV without name conflicts between packages. See, for example, the
#'   `kfold()` methods in the **rstanarm** and **brms** packages.
#'
#'   The **Value** section below describes the objects that `kfold()`
#'   methods should return in order to be compatible with
#'   [loo_compare()] and the **loo** package print methods.
#'
#'
#' @name kfold-generic
#' @param x A fitted model object.
#' @param ... Arguments to pass to specific methods.
#'
#' @return For developers defining a `kfold()` method for a class
#'   `"foo"`, the `kfold.foo()` function should return a list with class
#'   `c("kfold", "loo")` with at least the following named elements:
#'    * `"estimates"`: A `1x2` matrix containing the ELPD estimate and its
#'      standard error. The matrix must have row name "`elpd_kfold`" and column
#'      names `"Estimate"` and `"SE"`.
#'    * `"pointwise"`: A `Nx1` matrix with column name `"elpd_kfold"` containing
#'      the pointwise contributions for each data point.
#'
#'   It is important for the object to have at least these classes and
#'   components so that it is compatible with other functions like
#'   [loo_compare()] and `print()` methods.
#'
NULL

#' @rdname kfold-generic
#' @export
kfold <- function(x, ...) {
  UseMethod("kfold")
}

#' @rdname kfold-generic
#' @export
is.kfold <- function(x) {
  inherits(x, "kfold") && is.loo(x)
}

#' @export
dim.kfold <- function(x) {
  attr(x, "dims")
}
