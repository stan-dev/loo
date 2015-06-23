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
#' \item{\code{wcp = 20}}{the percentage of samples used for the generalized
#' Pareto fit estimate}
#' \item{\code{wtrunc = 3/4}}{for truncating very large weights to
#' \eqn{N}^\code{wtrunc} (set to zero for no truncation)}
#' \item{\code{fix_value = 100}}{if the largest value in any column of \code{lw
#' = -log_lik} is greater than the second largest by at least \code{fix_value},
#' then the second largest value will be set equal to the largest (set to zero
#' to skip this step)}
#'\item{\code{cores = \link[parallel]{detectCores}()}}{the number of cores to
#'      use for parallelization.}
#'}
#'
#' @return a named list. Returned for both LOO and WAIC are the expected log
#'   pointwise predictive density (\code{elpd} ), the estimated effective number
#'   of parameters (\code{p}), and the information criteria on the deviance scale
#'   (e.g. \code{looic = -2*elpd_loo}). Also returned are the pointwise
#'   contributions of each of these measures, standard errors, and the estimated
#'   shape parameter \eqn{k} for the Pareto fit to the importance ratios for
#'   each leave-one-out distribution.
#'
#' @seealso \code{\link{loo_and_waic_diff}}, \code{\link{loo-package}},
#'   \code{link{vgislw}}
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
    stop("'log_lik' should be a matrix")
  S <- nrow(log_lik)
  N <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  loo <- vgisloo(log_lik, ...)
  elpd_loo <- loo$loos
  p_loo <- lpd - elpd_loo
  looic <- -2 * elpd_loo
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  nms <- names(pointwise <- nlist(elpd_loo, p_loo, elpd_waic, p_waic, looic, waic))
  total <- unlist_lapply(pointwise, sum)
  se <- sqrt(N * unlist_lapply(pointwise, var))
  output <- as.list(c(total, se))
  names(output) <- c(nms, paste0("se_", nms))
  output$pointwise <- do.call("cbind", pointwise)
  output$pareto_k <- loo$ks
  output$info <- list(log_lik_nsims = S, log_lik_nobs = N)
  class(output) <- "loo"
  output
}
