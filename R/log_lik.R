#' Convenience function for extracting log-likelihood from
#' a \code{stanfit} object.
#'
#' @export
#' @param stanfit a \code{stanfit} (\pkg{rstan}) object.
#' @param parameter_name a character string naming the parameter (generated
#' quantity) in the Stan model corresponding to the log-likelihood.
#' @return a matrix of (post-warmup) extracted draws.
#' @seealso \code{\link[rstan]{stanfit-class}}
log_lik <- function(stanfit, parameter_name = "log_lik") {
  rstan_ok <- requireNamespace("rstan", quietly = TRUE)
  if (!rstan_ok) {
    stop("Please install the rstan package to use this function.")
  }
  rstan::extract(stanfit, parameter_name)[[parameter_name]]
}
