#' VGIS-LOO
#'
#' Very good importance sampling, log predictive densities
#'
#' @export
#'
#' @param log_lik an \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample and \eqn{N} is the number of data points (see
#'   \code{\link{extract_log_lik}}).
#' @param wcp the percentage of samples used for the genearlized Pareto fit
#'   estimate.
#' @param wtrunc for truncating very large weights to \eqn{N}^\code{wtrunc}. No
#'   trunction if \code{wtrunc=0}.
#' @param cores number of cores to use for parallelization (passed to
#'   \code{\link{vgislw}}).
#'
#' @note These functions are primarily intended for internal use, but are
#'   exported so that users can call them directly if desired. Users simply
#'   wishing to compute LOO and WAIC should use the \code{\link{loo_and_waic}}
#'   function.
#'
#' @return \code{vgisloo} calls \code{\link{vgislw}} and post-processes the
#'   results (a list with modified log weights and tail indices), returning a
#'   named list with components
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
vgisloo <- function(log_lik, wcp = 20, wtrunc = 3/4,
                    cores = parallel::detectCores()) {
  lw <- -1 * log_lik
  vgis <- vgislw(lw, wcp, wtrunc, cores)
  loos <- matrixStats::colLogSumExps(log_lik + vgis$lw)
  loo <- sum(loos)
  nlist(loo, loos, ks = vgis$k)
}
