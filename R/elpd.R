#' Generic (expected) log-predictive density (deprecated)
#'
#' As of loo 3.0.0, `elpd()` is **deprecated**. Please use [measure_elpd()]
#' instead. For full predictive performance workflows, see
#' [insample_pred_measure()] and [loo_pred_measure()].
#' See `vignette("migration-guide", package = "loo")` for a full mapping table.
#'
#' @details
#' The return type differs: `elpd()` returns class `"elpd_generic"` with
#' `elpd` and `ic` in `pointwise`; `measure_elpd()` returns class `"measure"`.
#'
#' The `elpd()` methods for arrays and matrices can compute the expected log
#' pointwise predictive density for a new dataset or the log pointwise
#' predictive density of the observed data (an overestimate of the elpd).
#' The `elpd()` function is an S3 generic and methods are provided for
#' 3-D pointwise log-likelihood arrays and matrices.
#'
#' @export
#' @param x A log-likelihood array or matrix. The **Methods (by class)**
#'   section, below, has detailed descriptions of how to specify the inputs for
#'   each method.
#' @param ... Currently ignored.
#'
#' @seealso [measure_elpd()], [insample_pred_measure()], [loo_pred_measure()],
#'   and the vignette *Holdout validation and K-fold cross-validation of Stan
#'   programs with the loo package*.
#'
#' @examples
#' \dontrun{
#' # Deprecated:
#' LLarr <- example_loglik_array()
#' elpd(LLarr)
#' # ->
#' measure_elpd(LLarr)
#' }
#'
elpd <- function(x, ...) {
  UseMethod("elpd")
}

#' @export
#' @templateVar fn elpd
#' @template array
#'
elpd.array <- function(x, ...) {
  .Deprecated("measure_elpd")
  ll <- llarray_to_matrix(x)
  .elpd_matrix_impl(ll)
}

#' @export
#' @templateVar fn elpd
#' @template matrix
#'
elpd.matrix <- function(x, ...) {
  .Deprecated("measure_elpd")
  .elpd_matrix_impl(x)
}


# internal ----------------------------------------------------------------
# used to avoid duplicated deprecation warning messages
.elpd_matrix_impl <- function(x) {
  pointwise <- pointwise_elpd_calcs(x)
  elpd_object(pointwise, dim(x))
}

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
