#' Generic (expected) log-predictive density
#'
#' The `elpd()` methods for arrays and matrices can compute the expected log
#' pointwise predictive density for a new dataset or the log pointwise
#' predictive density of the observed data (an overestimate of the elpd).
#'
#' @export
#' @param x A log-likelihood array or matrix. The **Methods (by class)**
#'   section, below, has detailed descriptions of how to specify the inputs for
#'   each method.
#' @param ... Currently ignored.
#'
#' @details The `elpd()` function is an S3 generic and methods are provided for
#'   3-D pointwise log-likelihood arrays and matrices.
#'
#' @seealso The vignette *Holdout validation and K-fold cross-validation of Stan
#'   programs with the loo package* for demonstrations of using the `elpd()`
#'   methods.
#'
#' @examples
#' # Calculate the lpd of the observed data
#' LLarr <- example_loglik_array()
#' elpd(LLarr)
#'
elpd <- function(x, ...) {
  UseMethod("elpd")
}

#' @export
#' @templateVar fn elpd
#' @template array
#'
elpd.array <- function(x, ...) {
  ll <- llarray_to_matrix(x)
  elpd.matrix(ll)
}

#' @export
#' @templateVar fn elpd
#' @template matrix
#'
elpd.matrix <- function(x, ...) {
  pointwise <- pointwise_elpd_calcs(x)
  elpd_object(pointwise, dim(x))
}



# internal ----------------------------------------------------------------
pointwise_elpd_calcs <- function(ll){
  elpd <- colLogSumExps(ll) - log(nrow(ll))
  ic <- -2 * elpd
  cbind(elpd, ic)
}

elpd_object <- function(pointwise, dims) {
  if (!is.matrix(pointwise)) stop("Internal error ('pointwise' must be a matrix)")

  cols_to_summarize <- colnames(pointwise)
  estimates <- table_of_estimates(pointwise[, cols_to_summarize, drop=FALSE])
  out <- nlist(estimates, pointwise)
  structure(
    out,
    dims = dims,
    class = c("elpd_generic", "loo")
  )
}
#' @export
print_dims.elpd_generic <- function(x, ...) {
  cat(
    "Computed from",
    paste(dim(x), collapse = " by "),
    "log-likelihood matrix using the generic elpd function\n"
  )
}
