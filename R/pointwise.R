#' Convenience function for extracting pointwise estimates
#'
#' @export
#' @param x A `loo` object, for example one returned by [loo()],
#'   [loo_subsample()], [loo_approximate_posterior()], [loo_moment_match()], etc.
#' @param estimate Which pointwise estimate to return. By default all are
#'   returned. The objects returned by the different functions ([loo()],
#'   [loo_subsample()], etc.) have slightly different estimates available.
#'   Typically at a minimum the estimates `elpd_loo`, `looic`, `mcse_elpd_loo`,
#'   `p_loo`, and `influence_pareto_k` will be available but there may be
#'   others.
#' @param ... Currently ignored.
#' @return If `estimate` is `NULL` then all pointwise estimates are returned in
#'   a matrix with one column per estimate and one row per observation.
#'   Otherwise a vector of length equal to the number of observations is
#'   returned containing the pointwise values for `estimate`.
#'
#' @examples
#' x <- loo(example_loglik_array())
#' head(pointwise(x))
#' pointwise(x, "elpd_loo")
#'
pointwise <- function(x, estimate = NULL, ...) {
  UseMethod("pointwise")
}

#' @rdname pointwise
#' @export
pointwise.loo <- function(x, estimate = NULL, ...) {
  pw <- x$pointwise
  if (is.null(pw)) {
    stop("No pointwise estimates found.", call. = FALSE)
  }
  if (is.null(estimate)) {
    return(pw)
  }
  estimates <- colnames(pw)
  if (!(estimate %in% estimates)) {
    stop(
      shQuote(estimate), " not found.",
      " Available estimates are: \n",
      paste(shQuote(estimates), collapse=", ")
    )
  }
  pw[, estimate]
}
