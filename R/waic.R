#' Widely applicable information criterion (WAIC)
#'
#' @export waic waic.matrix waic.function
#' @inheritParams loo
#' @param ... Other arguments. Currently ignored.
#'
#' @return A named list (of class \code{'loo'}) with components:
#'
#' \describe{
#' \item{\code{elpd_waic, se_elpd_waic}}{expected log pointwise predictive
#' density and standard error}
#' \item{\code{p_waic, se_p_waic}}{estimated effective number of parameters and
#' standard error}
#' \item{\code{waic, se_waic}}{\code{-2 * elpd_waic} (i.e., converted to the
#' deviance scale) and standard error}
#' \item{\code{pointwise}}{the pointwise contributions of each of the above
#' measures}
#' }
#'
#' @seealso \code{\link{compare}}, \code{\link{print.loo}},
#' \code{\link{loo-package}}
#'
#'
#' @examples
#' \dontrun{
#' log_lik1 <- extract_log_lik(stanfit1)
#' waic1 <- waic(log_lik1)
#' waic1
#'
#' log_lik2 <- extract_log_lik(stanfit2)
#' waic2 <- waic(log_lik2)
#' print(waic2, digits = 4)
#'
#' waic_diff <- compare(waic1, waic2)
#' print(waic_diff, digits = 2)
#' }
#'
waic <- function(x, ...) {
  UseMethod("waic")
}

#' @describeIn waic
#'
#' An \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the posterior
#' sample (the number of simulations) and \eqn{N} is the number of data points.
#' Typically (but not restricted to be) the object returned by
#' \code{\link{extract_log_lik}}.
#'
#' @export
#'
waic.matrix <- function(x, ...) {
  out <- pointwise_waic(log_lik = x)
  structure(out, log_lik_dim = dim(x), class = "loo")
}

#' @describeIn waic
#'
#'  A function \code{f} that takes arguments \code{i}, \code{data},
#'  and \code{draws} and returns a vector containing the log-likelihood for
#'  the \code{i}th observation evaluated at each posterior draw.
#'
#'  The \code{args} argument must also be specified and should be a named list
#'  with the following components:
#'  \itemize{
#'    \item \code{draws}: An object containing the posterior draws for any
#'    parameters needed to compute the pointwise log-likelihood.
#'    \item \code{data}: A matrix of data (e.g. observed outcome and predictors)
#'    needed to compute the pointwise log-likelihood and in the appropriate form
#'    so that \code{f(i=i, data=data[i,,drop=FALSE], draws=draws)} returns the
#'    vector of log-likelihods for the \code{i}th observation.
#'    \item \code{N}: The number of observations.
#'    \item \code{S}: The size of the posterior sample.
#'  }
#'
#' @export
#'
waic.function <- function(x, ..., args) {
  out <- pointwise_waic(llfun = x, llargs = args)
  structure(out, log_lik_dim = with(args, c(S,N)), class = "loo")
}
