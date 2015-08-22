#' @describeIn <%= fn %>
#'
#'  A function \eqn{f} that takes arguments \code{i}, \code{data},
#'  and \code{draws} and returns a vector containing the log-likelihood for
#'  the \code{i}th observation evaluated at each posterior draw.
#'
#'  The \code{args} argument must also be specified and should be a named list
#'  with the following components:
#'  \itemize{
#'    \item \code{draws}: An object containing the posterior draws for any
#'    parameters needed to compute the pointwise log-likelihood.
#'    \item \code{data}: An object containing data (e.g. observed outcome and
#'    predictors) needed to compute the pointwise log-likelihood. \code{data}
#'    should be in an appropriate form so that \eqn{f}\code{(i=i,
#'    data=data[i,,drop=FALSE], draws=draws)} returns the \code{S}-vector of
#'    log-likelihoods for the \code{i}th observation.
#'    \item \code{N}: The number of observations.
#'    \item \code{S}: The size of the posterior sample.
#'  }
#'
