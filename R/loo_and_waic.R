#' Approximate LOO-CV and WAIC for Bayesian models
#'
#' @export
#' @param log_lik an \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample (the number of simulations) and \eqn{N} is the number of
#'   data points. Typically (but not restricted to be) the object returned by
#'   \code{\link{extract_log_lik}}.
#' @param ... optional arguments to pass to \code{\link{vgislw}}. Possible
#' arguments and their defaults are:
#' \describe{
#' \item{\code{wcp = 0.2}}{the proportion of importance weights to use for the
#' generalized Pareto fit. The \code{100*wcp}\% largest weights are used as the
#' sample from which to estimate the parameters \eqn{k} and \eqn{\sigma} of
#' the generalized Pareto distribution.}
#' \item{\code{thresh = 100}}{when used for the generalized Pareto fit, the
#' largest \code{100*wcp}\% of the importance weights are modified according to
#' \code{pmax(x, max(x) - thresh)} for numerical stability.}
#' \item{\code{kmax = 2}}{the maximum allowed value for the generalized Pareto
#' shape parameter \eqn{k}.}
#' \item{\code{wtrunc = 3/4}}{for truncating very large weights to
#' \eqn{S}^\code{wtrunc} (set to zero for no truncation).}
#'\item{\code{cores = \link[parallel]{detectCores}()}}{the number of cores to
#'      use for parallelization.}
#'}
#'
#' We recommend using the default values for the \code{vgislw} arguments unless
#' there are problems (e.g. \code{NA} or \code{NaN} results).
#'
#' @return a named list with class \code{'loo'}.
#'
#' Returned for both LOO and WAIC are the expected log pointwise predictive
#' density (\code{elpd}), the estimated effective number of parameters
#' (\code{p}), and the information criteria on the deviance scale (e.g.
#' \code{looic = -2*elpd_loo}). Also returned are the pointwise contributions of
#' each of these measures, standard errors, and the estimated shape parameter
#' \eqn{k} for the Pareto fit to the importance ratios for each leave-one-out
#' distribution.
#'
#' @seealso \code{\link{loo-package}}, \code{\link{print.loo}},
#' \code{\link{loo_and_waic_diff}}
#'
#' @examples
#' \dontrun{
#' log_lik <- extract_log_lik(stanfit)
#' loo <- loo_and_waic(log_lik)
#' print(loo, digits = 3)
#' }
#'
#' @importFrom matrixStats colVars
#'
loo_and_waic <- function(log_lik, ...) {
  if (!is.matrix(log_lik))
    stop('log_lik should be a matrix')
  loo <- vgisloo(log_lik, ...)
  lpd <- logColMeansExp(log_lik)
  elpd_loo <- loo$loos
  p_loo <- lpd - elpd_loo
  looic <- -2 * elpd_loo
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  nms <- names(pointwise <- nlist(elpd_loo, p_loo, elpd_waic, p_waic, looic, waic))
  total <- unlist_lapply(pointwise, sum)
  se <- sqrt(ncol(log_lik) * unlist_lapply(pointwise, var))
  output <- as.list(c(total, se))
  names(output) <- c(nms, paste0("se_", nms))
  output$pointwise <- cbind_list(pointwise)
  output$pareto_k <- loo$pareto_k
  attr(output, "log_lik_dim") <- dim(log_lik)
  class(output) <- "loo"
  output
}
