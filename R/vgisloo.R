#' Very good importance sampling LOO log predictive densities
#'
#' @export
#' @param log_lik an \eqn{s} by \eqn{n} matrix, where \eqn{s} is the size of the
#'   posterior sample and \eqn{n} is the number of data points (see
#'   \code{\link{extract_log_lik}}).
#' @param wcp the percentage of samples used for the genearlized Pareto fit
#'   estimate
#' @param wtrunc for truncating very large weights to \eqn{n}^\code{wtrunc}. No
#'   trunction if \code{wtrunc=0}.
#' @param cores number of cores to use for parallelization.
#'
#' @return A list with components \describe{ \item{\code{loo}}{the sum of the
#'   LOO log predictive densities} \item{\code{loos}}{the individual LOO log
#'   predictive densities} \item{\code{ks}}{the estimate of the tail indices} }
#'

vgisloo <- function(log_lik, wcp = 20, wtrunc = 3/4, cores = parallel::detectCores()) {
  lw <- -1 * log_lik
  vgis <- vgislw(lw, wcp, wtrunc, cores)
  loos <- matrixStats::colLogSumExps(log_lik + vgis$lw)
  loo <- sum(loos)
  nlist(loo, loos, ks = vgis$k)
}
