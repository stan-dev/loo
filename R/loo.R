#' Leave-one-out cross-validation (LOO)
#'
#' Efficient approximate leave-one-out cross-validation
#'
#' @export loo loo.matrix loo.function
#' @param x A log-likelihood matrix or function. See the \strong{Methods (by
#'   class)} section below for a detailed description.
#' @param args Only required if \code{x} is a function. A list containing
#'   the data required to specify the arguments to the function. See the
#'   \strong{Methods (by class)} section below for how \code{args} should be
#'   specified.
#' @param ... Optional arguments to pass to \code{\link{psislw}}. Possible
#'   arguments and their defaults are:
#' \describe{
#'  \item{\code{wcp = 0.2}}{The proportion of importance weights to use for the
#'    generalized Pareto fit. The \code{100*wcp}\% largest weights are used as the
#'    sample from which to estimate the parameters \eqn{k} and \eqn{\sigma} of
#'    the generalized Pareto distribution.}
#'  \item{\code{wtrunc = 3/4}}{For truncating very large weights to
#'    \eqn{S}^\code{wtrunc} (set to zero for no truncation).}
#'  \item{\code{cores = \link[parallel]{detectCores}()}}{The number of cores to
#'    use for parallelization.}
#'}
#'
#'  We recommend using the default values for the \code{psislw} arguments unless
#'  there are problems (e.g. \code{NA} or \code{NaN} results).
#'
#' @return A named list with class \code{'loo'} and components:
#'
#' \describe{
#'  \item{\code{elpd_loo, se_elpd_loo}}{expected log pointwise predictive density
#'    and standard error}
#'  \item{\code{p_loo, se_p_loo}}{estimated effective number of parameters and
#'    standard error}
#'  \item{\code{looic, se_looic}}{\code{-2 * elpd_loo} (i.e., converted to the
#'    deviance scale) and standard error}
#'  \item{\code{pointwise}}{a matrix containing the pointwise contributions of each
#'    of the above measures}
#'  \item{\code{pareto_k}}{a vector containing the estimates of the shape
#'    parameter \eqn{k} for the generaelized Pareto fit to the importance ratios
#'    for each leave-one-out distribution. See PSIS-LOO section in
#'    \code{\link{loo-package}} for details about interpreting \eqn{k}. (Also, by
#'    default, the print method for \code{'loo'} objects will provide warnings
#'    about problematic values of \eqn{k}. It is also possible to plot the values
#'    of \eqn{k} by setting the optional argument \code{plot_k} to \code{TRUE} (the
#'    default is not to plot). See \code{\link{print.loo}}.)}
#' }
#'
#' @seealso \code{\link{loo-package}}, \code{\link{print.loo}},
#' \code{\link{compare}}
#'
#' @examples
#' \dontrun{
#' log_lik1 <- extract_log_lik(stanfit1)
#' loo1 <- loo(log_lik1)
#' print(loo1, digits = 3)
#'
#' log_lik2 <- extract_log_lik(stanfit2)
#' loo2 <- loo(log_lik2)
#' loo2
#'
#' loo_diff <- compare(loo1, loo2)
#' loo_diff
#' print(loo_diff, digits = 5)
#'}
#'
loo <- function(x, ...) {
  UseMethod("loo")
}

#' @describeIn loo
#'
#' An \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the posterior
#' sample (the number of simulations) and \eqn{N} is the number of data points.
#' Typically (but not restricted to be) the object returned by
#' \code{\link{extract_log_lik}}.
#'
#' @export
#'
loo.matrix <- function(x, ...) {
  if (any(is.na(x)))
    stop("NA log-likelihood values found.", call. = FALSE)
  psis <- psislw(lw = -1 * x, ...)
  out <- pointwise_loo(psis, x)
  structure(out, log_lik_dim = dim(x), class = "loo")
}

#' @describeIn loo
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
#'
#'
#'
#' @export
#'
loo.function <- function(x, ..., args) {
  if (missing(args)) stop("args must be specified", call. = FALSE)
  psis <- psislw(..., llfun = x, llargs = args)
  out <- pointwise_loo(psis = psis, llfun = x, llargs = args)
  structure(out, log_lik_dim = with(args, c(S,N)), class = "loo")
}

