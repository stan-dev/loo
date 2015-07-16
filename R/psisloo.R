#' psisloo
#'
#' \code{\link{psislw}} performs Pareto smooth importance sampling. \code{psisloo}
#' calls \code{psislw} and does useful post-processing for computing LOO and
#' WAIC.
#'
#' @keywords internal
#' @export
#'
#' @param log_lik an \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample and \eqn{N} is the number of data points (see
#'   \code{\link{extract_log_lik}}).
#' @param ... optional arguments to pass to \code{\link{psislw}}.
#'
#' @note This function is primarily intended for internal use, but is
#'   exported so that users can call it directly if desired. Users simply
#'   wishing to compute LOO should use the \code{\link{loo}} function.
#'
#' @return \code{psislw} returns a list with modified log weights and tail
#'   indices. \code{psisloo} calls \code{psislw(lw = -log_lik,...)},
#'   post-processes the results, and returns a named list with components
#' \describe{
#' \item{\code{loos}}{the individual LOO log predictive densities}
#' \item{\code{pareto_k}}{shape parameter estimates for the generalized
#' Pareto distribution.}
#' }
#'
#' @details See the 'PSIS-LOO' section in \code{\link{loo-package}}.
#'
#' @seealso \code{\link{psislw}}, \code{\link{loo}},
#' \code{\link{loo-package}}.
#'
#' @importFrom matrixStats colLogSumExps
#'
psisloo <- function(log_lik, ...) {
  lw <- -1 * log_lik
  psis <- psislw(lw, ...)
  loos <- colLogSumExps(log_lik + psis$lw_smooth)
  nlist(loos, pareto_k = psis$pareto_k)
}
