#' Approximate LOO-CV and WAIC for Bayesian models
#'
#' @export
#' @param log_lik an \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample (the number of simulations) and \eqn{N} is the number of
#'   data points. Typically (but not restricted to be) the object returned by
#'   \code{\link{extract_log_lik}}.
#' @param ... optional arguments to pass to \code{\link{vgislw}} (via
#' \code{\link{vgisloo}}). See \code{\link{vgislw}}.
#' @return a named list. Returned for both LOO and WAIC are the expected log
#'   pointwise predictive density (\code{elpd} ), the estimated effective number
#'   of parameters (\code{p}), and the information criteria on the deviance scale
#'   (e.g. \code{looic = -2*elpd_loo}). Also returned are the pointwise
#'   contributions of each of these measures, standard errors, and the estimated
#'   shape parameter \eqn{k} for the Pareto fit to the importance ratios for
#'   each leave-one-out distribution.
#'
#' @details The possible arguments that can be specified in \code{...} are
#' wcp, wtrunc, fix_value, cores
#'
#' @seealso \code{\link{loo_and_waic_diff}}, \code{\link{loo-package}}
#' @examples
#' \dontrun{
#' log_lik <- extract_log_lik(stanfit)
#' loo <- loo_and_waic(log_lik)
#' print(loo, digits = 3)
#' }
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
  p_waic <- matrixStats::colVars(log_lik)
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
