#' @describeIn <%= fn %>
#'   A function \code{f} that takes arguments \code{data_i} and \code{draws} and
#'   returns a vector containing the log-likelihood for a single observation
#'   \code{i} evaluated at each posterior draw. The function should be written
#'   such that, for each observation \code{i} in \code{1:N}, evaluating
#'   \code{f(data_i = data[i,, drop=FALSE], draws = draws)} results in a vector
#'   of length \code{S} (size of posterior sample).
#'
#'  If using the function method then the arguments \code{data}
#'  and \code{draws} must also be specified in the call to \code{loo}:
#'  \itemize{
#'    \item \code{data}: An object containing data (e.g. observed outcome and
#'    predictors) needed to compute the pointwise log-likelihood. \code{data}
#'    must be an object that can be indexed as described above
#'    (\code{data[i,,drop=FALSE]}) to access the data for the \code{i}th
#'    observation. \code{data} must have a \code{dim} attribute such that
#'    \code{dim(data)[1]} is equal to \code{N}, the number of observations.
#'    \item \code{draws}: An object containing the posterior draws for any
#'    parameters needed to compute the pointwise log-likelihood. \code{draws}
#'    must have a \code{dim} attribute such that \code{dim(draws)[1]} is equal
#'    to \code{S}, the size of the posterior sample.
#'  }
#'
