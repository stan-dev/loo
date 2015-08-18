#' Extract log-likelihood from a Stan model
#'
#' Convenience function for extracting the pointwise log-likelihood from a
#' fitted Stan model.
#'
#' @export
#' @param stanfit A \code{stanfit} object (\pkg{rstan} package).
#' @param parameter_name A character string naming the parameter (or generated
#'   quantity) in the Stan model corresponding to the log-likelihood.
#' @return An \eqn{S} by \eqn{N} matrix of (post-warmup) extracted draws, where
#'   \eqn{S} is the number of simulations and \eqn{N} is the number of data
#'   points.
#'
#' @details Stan does not automatically compute and store the log-likelihood. It
#'   is up to the user to incorporate it into the Stan program if it is to be
#'   extracted after fitting the model. In a Stan model, the pointwise log
#'   likelihood can be coded as a vector in the transformed parameters block
#'   (and then summed up in the model block) or it can be coded entirely in the
#'   generated quantities block. We recommend using the generated quantities
#'   block so that the computations are carried out only once per iteration
#'   rather than once per HMC leapfrog step.
#'
#'   For example, the following is the \code{generated quantities} block for
#'   computing and saving the log-likelihood for a linear regression model with
#'   \code{N} data points, outcome \code{y}, predictor matrix \code{X},
#'   coefficients \code{beta}, and standard deviation \code{sigma}:
#'
#'  \code{vector[N] log_lik;}
#'
#'  \code{for (n in 1:N) log_lik[n] <- normal_log(y[n], X[n] * beta, sigma);}
#'
#' @references
#' Stan Development Team (2015). Stan: A C++ library for probability and
#' sampling, version 2.6. \url{mc-stan.org}.
#'
#' Stan Development Team (2015). RStan, version 2.6.
#' \url{mc-stan.org/rstan.html}.
#'
#' @examples
#' \dontrun{
#' log_lik <- extract_log_lik(stanfit, "log_lik")
#' }
#'
extract_log_lik <- function(stanfit, parameter_name = "log_lik") {
  if (!inherits(stanfit, "stanfit"))
    stop("Not a stanfit object.")
  if (stanfit@mode != 0)
    stop("Stan model does not contain samples.")
  posterior <- as.matrix(stanfit)
  nms <- colnames(posterior)
  pattern <- paste0("^", parameter_name, "\\[")
  keep <- grep(pattern, nms)
  log_lik <- posterior[, keep, drop = FALSE]
  colnames(log_lik) <- NULL
  log_lik
}

