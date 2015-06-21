#' Convenience function for extracting the point-wise log-likelihood from a
#' \code{stanfit} object
#'
#' @export
#' @param stanfit a \code{stanfit} (\pkg{rstan}) object.
#' @param parameter_name a character string naming the parameter (generated
#'   quantity) in the Stan model corresponding to the log-likelihood.
#' @return a matrix of (post-warmup) extracted draws.
#'
#' @details In the Stan program the point-wise log likelihood can either be
#'   defined as a vector in the transformed parameters block (and then summed up
#'   in the model block) or it can be coded entirely in the generated quantities
#'   block using Stan's log probability density functions. The result after
#'   running Stan will be an \eqn{S} by \eqn{N} matrix of point-wise
#'   log-likelihood values, where \eqn{S} is the size of the posterior sample
#'   (the number of draws/simulations) and \eqn{N} is the number of data points.
#'
#' @note The \pkg{rstan} package is required in order to use this function
#' @seealso \code{\link[rstan]{stanfit-class}}
#' @examples
#' \dontrun{
#' log_lik <- extract_log_lik(stanfit, "log_lik")
#' }
#'
extract_log_lik <- function(stanfit, parameter_name = "log_lik") {
  rstan_ok <- requireNamespace("rstan", quietly = TRUE)
  if (!rstan_ok) {
    stop("Please install the rstan package to use this function.")
  }
  rstan::extract(stanfit, parameter_name)[[parameter_name]]
}
