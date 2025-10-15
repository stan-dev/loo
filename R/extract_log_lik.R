#' Extract pointwise log-likelihood from a Stan model
#'
#' Convenience function for extracting the pointwise log-likelihood
#' matrix or array from a `stanfit` object from the \pkg{rstan} package.
#' Note: recent versions of \pkg{rstan} now include a `loo()` method for
#' `stanfit` objects that handles this internally.
#'
#' @export
#' @param stanfit A `stanfit` object (\pkg{rstan} package).
#' @param parameter_name A character string naming the parameter (or generated
#'   quantity) in the Stan model corresponding to the log-likelihood.
#' @param merge_chains If `TRUE` (the default), all Markov chains are
#'   merged together (i.e., stacked) and a matrix is returned. If `FALSE`
#'   they are kept separate and an array is returned.
#' @return If `merge_chains=TRUE`, an \eqn{S} by \eqn{N} matrix of
#'   (post-warmup) extracted draws, where \eqn{S} is the size of the posterior
#'   sample and \eqn{N} is the number of data points. If
#'   `merge_chains=FALSE`, an \eqn{I} by \eqn{C} by \eqn{N} array, where
#'   \eqn{I \times C = S}{I * C = S}.
#'
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
#'   For example, the following is the `generated quantities` block for
#'   computing and saving the log-likelihood for a linear regression model with
#'   `N` data points, outcome `y`, predictor matrix `X`,
#'   coefficients `beta`, and standard deviation `sigma`:
#'
#'  `vector[N] log_lik;`
#'
#'  `for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | X[n, ] * beta, sigma);`
#'
#' @references
#' Stan Development Team (2017). The Stan C++ Library, Version 2.16.0.
#' <https://mc-stan.org/>
#'
#' Stan Development Team (2017). RStan: the R interface to Stan, Version 2.16.1.
#' <https://mc-stan.org/>
#'
extract_log_lik <-
  function(stanfit,
           parameter_name = "log_lik",
           merge_chains = TRUE) {
    if (!inherits(stanfit, "stanfit"))
      stop("Not a stanfit object.", call. = FALSE)
    if (stanfit@mode != 0)
      stop("Stan model does not contain posterior draws.", call. = FALSE)
    if (!requireNamespace("rstan", quietly = TRUE))
      stop("Please load the 'rstan' package.", call. = FALSE)

    if (merge_chains) {
      log_lik <- as.matrix(stanfit, pars = parameter_name)
    } else {
      log_lik <- as.array(stanfit, pars = parameter_name)
    }
    unname(log_lik)
  }
