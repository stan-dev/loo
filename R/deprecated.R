#' Approximate LOO-CV and WAIC for Bayesian models
#'
#' This function is deprecated. Please use \code{\link{loo}} or
#' \code{\link{waic}} instead.
#'
#' @keywords internal
#' @export
#'
#' @param log_lik an \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample (the number of simulations) and \eqn{N} is the number of
#'   data points. Typically (but not restricted to be) the object returned by
#'   \code{\link{extract_log_lik}}.
#' @param ... optional arguments to pass to \code{\link{psislw}}. Possible
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
#' We recommend using the default values for the \code{psislw} arguments unless
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
#' @importFrom matrixStats colVars
#'
loo_and_waic <- function(log_lik, ...) {

  .Deprecated("loo() or waic()")

  if (!is.matrix(log_lik))
    stop('log_lik should be a matrix')
  loo <- psisloo(log_lik, ...)
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


#' Model comparison
#'
#' This function is deprecated. Please use \code{\link{compare}} instead.
#'
#' @keywords internal
#' @export
#'
#' @param loo1,loo2 lists (of class \code{'loo'}).
#' @return a named list with class \code{'compare.loo'}.
#'
#' @details When comparing two fitted models, we can estimate the difference in
#'   their expected predictive accuracy by the difference in \code{elpd_waic} or
#'   \code{elpd_loo} (multiplied by -2, if desired, to be on the deviance
#'   scale). To compute the standard error of this difference we can use a
#'   paired estimate to take advantage of the fact that the same set of \eqn{N}
#'   data points is being used to fit both models. We would think that these
#'   calculations would be most useful when \eqn{N} is large, because then
#'   non-normality of the distribution is not such an issue when estimating the
#'   uncertainty of these sums. In any case, we suspect that these standard
#'   errors, for all their flaws, should give a better sense of uncertainty than
#'   what is obtained using the current standard approach of comparing
#'   differences of deviances to a Chi-squared distribution, a practice derived
#'   for Gaussian linear models or asymptotically and which only applies to
#'   nested models in any case.
#'
#'
loo_and_waic_diff <- function(loo1, loo2) {

  .Deprecated("compare()")

  p1 <- loo1$pointwise
  p2 <- loo2$pointwise
  N1 <- nrow(p1)
  N2 <- nrow(p2)
  if (N1 != N2) {
    msg <- paste("Models should have the same number of data points.",
                 "Found N1 =", N1, "and N2 =", N2)
    stop(msg)
  }
  sqrtN <- sqrt(N1)
  loo_diff <- p2[, "elpd_loo"] - p1[, "elpd_loo"]
  waic_diff <- p2[, "elpd_waic"] - p1[, "elpd_waic"]
  diff <- list(elpd_loo_diff = sum(loo_diff),
               lpd_loo_diff = sqrtN * sd(loo_diff),
               elpd_waic_diff = sum(waic_diff),
               lpd_waic_diff = sqrtN * sd(waic_diff))
  class(diff) <- "compare.loo"
  diff
}

