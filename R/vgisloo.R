#' vgisloo
#' @noRd
#'
#' \code{vgislw} performs very good importance sampling. \code{vgisloo} calls
#' \code{vgislw} and does useful post-processing for computing LOO and WAIC.
#'
#' @param log_lik an \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample and \eqn{N} is the number of data points (see
#'   \code{\link{extract_log_lik}}).
#' @param ... optional arguments to pass to \code{\link{vgislw}}.
#'
#' @note This function is primarily intended for internal use, but is
#'   exported so that users can call it directly if desired. Users simply
#'   wishing to compute LOO and WAIC should use the \code{\link{loo_and_waic}}
#'   function.
#'
#' @return \code{vgislw} returns a list with modified log weights and tail
#'   indices. \code{vgisloo} calls \code{vgislw(lw = -log_lik,...)},
#'   post-processes the results, and returns a named list with components
#' \describe{
#' \item{\code{loo}}{the sum of the LOO log predictive densities}
#' \item{\code{loos}}{the individual LOO log predictive densities}
#' \item{\code{ks}}{the estimate of the tail indices.}
#' }
#'
#' @seealso \code{\link{loo_and_waic}}, \code{\link{loo-package}}
#' @references
#' Vehtari, A., and Gelman, A. (2015). Very good importance sampling.
#'
#' @importFrom matrixStats colLogSumExps
#'
vgisloo <- function(log_lik, ...) {
  lw <- -1 * log_lik
  vgis <- vgislw(lw, ...)
  loos <- colLogSumExps(log_lik + vgis$lw)
  loo <- sum(loos)
  nlist(loo, loos, ks = vgis$k)
}
