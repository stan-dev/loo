#' Identify problematic observations
#'
#' Find observations for which the estimated Pareto shape parameter \eqn{k} is
#' larger than some \code{threshold} value. See the PSIS-LOO section in
#' \code{\link{loo-package}} for details about the interpretation of \eqn{k}.
#'
#' @export
#' @param x An object created by \code{\link{loo}}.
#' @param threshold The threshold value.
#' @return An integer vector indicating which observations have Pareto \eqn{k}
#'   estimates above \code{threshold}.
#'
pareto_k_ids <- function(x, threshold = 0.5) {
  if (is.null(x$pareto_k))
    stop("No Pareto k estimates found.", call. = FALSE)
  which(x$pareto_k > threshold)
}
