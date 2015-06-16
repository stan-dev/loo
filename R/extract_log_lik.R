#' Convenience function for extracting log-likelihood from a \code{stanfit}
#' object.
#'
#' @export
#' @param object a \code{stanfit} (\pkg{rstan}) object.
#' @param ... other arguments. Currently only \code{parameter_name}, which
#'   defaults to \code{'log_lik'}, and should be a character string naming the
#'   parameter (generated quantity) in the Stan model corresponding to the
#'   log-likelihood.
#' @return a matrix of (post-warmup) extracted draws.
#' @seealso \code{\link[rstan]{stanfit-class}}
#' @examples
#' \dontrun{
#' log_lik <- extract_log_lik(stanfit, "log_lik")
#' }
extract_log_lik <- function(stanfit, parameter_name = "log_lik") {
  rstan_ok <- requireNamespace("rstan", quietly = TRUE)
  if (!rstan_ok) {
    stop("Please install the rstan package to use this function.")
  }
  rstan::extract(stanfit, parameter_name)[[parameter_name]]
}
