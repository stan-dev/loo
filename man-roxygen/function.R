#' @describeIn <%= fn %>
#'   A function \code{f} that takes arguments \code{data_i} and \code{draws} and
#'   returns a vector containing the log-likelihood for a single observation
#'   \code{i} evaluated at each posterior draw. The function should be written
#'   such that, for each observation \code{i} in \code{1:N}, evaluating
#'   \code{f(data_i = data[i,, drop=FALSE], draws = draws)} results in a vector
#'   of length \code{S} (size of posterior sample). The log-likelihood function
#'   can also have additional arguments but \code{data_i} and \code{draws} are
#'   required.
#'
#'  If using the function method then the arguments \code{data}
#'  and \code{draws} must also be specified in the call to \code{loo}:
#'  \itemize{
#'    \item \code{data}: A data frame or matrix containing the data (e.g.
#'    observed outcome and predictors) needed to compute the pointwise
#'    log-likelihood. For each observation \code{i}, the \code{i}th row of
#'    \code{data} will be passed to the \code{data_i} argument of the
#'    log-likelihood function.
#'    \item \code{draws}: An object containing the posterior draws for any
#'    parameters needed to compute the pointwise log-likelihood. Unlike
#'    \code{data}, which is indexed by observation, for each observation the
#'    entire object \code{draws} will be passed to the \code{draws} argument of
#'    the log-likelihood function.
#'    \item The \code{...} can be used to pass additional arguments to your
#'    log-likelihood function. These arguments ase used like the \code{draws}
#'    argument in that they are recycled for each observation.
#'  }
#'
