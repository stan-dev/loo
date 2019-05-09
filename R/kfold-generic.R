#' Generic function for K-fold cross-validation for developers
#'
#' @description For developers of Bayesian modeling packages, \pkg{loo} includes
#'   a generic function \code{kfold} so that methods may be defined for K-fold
#'   CV without name conflicts between packages. See, for example, the
#'   \code{kfold} methods in \pkg{rstanarm} and \pkg{brms}.
#'
#'   The \strong{Value} section below describes the objects that \code{kfold}
#'   methods should return in order to be compatible with
#'   \code{\link{loo_compare}} and the \pkg{loo} package print methods.
#'
#'
#' @name kfold-generic
#' @param x A fitted model object.
#' @param ... Arguments to pass to specific methods.
#'
#' @return For developers defining a \code{kfold} method for a class
#'   \code{"foo"}, the \code{kfold.foo} function should return a list with class
#'   \code{c("kfold", "loo")} with at least the elements
#'   \itemize{
#'    \item \code{"estimates"}: a 1x2 matrix with column names "Estimate" and "SE"
#'    containing the ELPD estimate and its standard error.
#'    \item \code{"pointwise"}: an Nx1 matrix with column name "elpd_kfold" containing
#'    the pointwise contributions for each data point.
#'   }
#'   It is important for the object to have at least these classes and
#'   components so that it is compatible with other functions like
#'   \code{\link{loo_compare}} and print methods.
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
