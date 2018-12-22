#' Generic function for K-fold cross-validation for developers
#'
#' For developers of modeling packages, \pkg{loo} includes a generic function
#' \code{kfold} so that methods may be defined for K-fold CV without name
#' conflicts between packages. See, e.g., the \code{kfold.stanreg} method in
#' \pkg{rstanarm} and the \code{kfold.brmsfit} method in \pkg{brms}.
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
