#' LOO and WAIC
#'
#' @export
#' @param log_lik an nsims by nobs matrix, typically (but not restricted to be)
#' the object returned by \code{\link{log_lik}}.
#' @param cores number of cores to use for parallization.
#' @return a list of class \code{'loo'}.
#'
#' @details Leave-one-out cross-validation (LOO) and the widely applicable
#' information criterion (WAIC) are methods for estimating pointwise out-of-sample
#' prediction accuracy from a fitted Bayesian model using the log-likelihood
#' evaluated at the posterior simulations of the parameter values. LOO and WAIC
#' have various advantages over simpler estimates of predictive error such as
#' AIC and DIC but are less used in practice because they involve additional
#' computational steps. Here we lay out fast and stable computations for LOO and
#' WAIC that can be performed using existing simulation draws. We compute LOO
#' using very good importance sampling (VGIS), a new procedure for regularizing
#' importance weights. As a byproduct of our calculations, we also obtain
#' approximate standard errors for estimated predictive errors and for comparing
#' of pre- dictive errors between two models.
#'
loo_and_waic <- function(log_lik, cores = parallel::detectCores()) {
  # log_lik should be a matrix with nrow = nsims and ncol = nobs

  if (!is.matrix(log_lik))
    stop("'log_lik' should be a matrix")
  S <- nrow(log_lik)
  N <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  loo <- vgisloo(log_lik, cores)
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
