#' Leave-one-out cross-validation (LOO)
#'
#' Efficient approximate leave-one-out cross-validation for Bayesian models.
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
#'  \item{\code{cores = getOption("loo.cores", parallel::detectCores())}}{The
#'  number of cores to use for parallelization. This can be set for an entire R
#'  session by \code{options(loo.cores = NUMBER)}. The default is
#'  \code{\link[parallel]{detectCores}}().}
#'}
#'
#'  We recommend using the default values for the \code{psislw} arguments unless
#'  there are problems (e.g. \code{NA} or \code{NaN} results).
#'
#' @return A named list with class \code{'loo'} and components:
#'
#' \describe{
#'  \item{\code{elpd_loo, se_elpd_loo}}{Expected log pointwise predictive density
#'    and standard error.}
#'  \item{\code{p_loo, se_p_loo}}{Estimated effective number of parameters and
#'    standard error.}
#'  \item{\code{looic, se_looic}}{The LOO information criterion
#'    (\code{-2*elpd_loo}, i.e., converted to deviance scale) and standard
#'    error.}
#'  \item{\code{pointwise}}{A matrix containing the pointwise contributions of each
#'    of the above measures.}
#'  \item{\code{pareto_k}}{A vector containing the estimates of the shape
#'    parameter \eqn{k} for the generaelized Pareto fit to the importance ratios
#'    for each leave-one-out distribution. See PSIS-LOO section in
#'    \code{\link{loo-package}} for details about interpreting \eqn{k}.
#'    (By default, the \code{\link[=print.loo]{print}} method for \code{'loo'}
#'    objects will also provide warnings about problematic values of \eqn{k}.)}
#' }
#'
#' @seealso \code{\link{loo-package}}, \code{\link{print.loo}},
#' \code{\link{compare}}
#'
#' @examples
#' \dontrun{
#' ### Usage with stanfit objects
#' log_lik1 <- extract_log_lik(stanfit1) # see ?extract_log_lik
#' loo1 <- loo(log_lik1)
#' print(loo1, digits = 3)
#'
#' log_lik2 <- extract_log_lik(stanfit2)
#' (loo2 <- loo(log_lik2))
#' compare(loo1, loo2)
#' }
#'
#' ### Using log-likelihood function instead of matrix
#' set.seed(024)
#'
#' # Simulate data and draw from posterior
#' N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
#' p <- rbeta(1, a0, b0)
#' y <- rbinom(N, size = K, prob = p)
#' a <- a0 + sum(y); b <- b0 + N * K - sum(y)
#' draws <- rbeta(S, a, b)
#' data <- data.frame(y,K)
#'
#' llfun <- function(i, data, draws) {
#'   dbinom(data$y, size = data$K, prob = draws, log = TRUE)
#' }
#' loo_with_fn <- loo(llfun, args = nlist(data, draws, N, S), cores = 1)
#'
#' # Check that we get same answer if using log-likelihood matrix
#' log_lik_mat <- sapply(1:N, function(i) llfun(i, data[i,, drop=FALSE], draws))
#' loo_with_mat <- loo(log_lik_mat, cores = 1)
#' all.equal(loo_with_mat, loo_with_fn)
#'
loo <- function(x, ...) {
  UseMethod("loo")
}

#' @export
#' @templateVar fn loo
#' @template matrix
#'
loo.matrix <- function(x, ...) {
  if (any(is.na(x))) stop("NA log-likelihood values found.")
  psis <- psislw(lw = -1 * x, ..., COMPUTE_LOOS = TRUE)
  out <- pointwise_loo(psis, x)
  structure(out, log_lik_dim = dim(x), class = "loo")
}

#' @export
#' @templateVar fn loo
#' @template function
#'
loo.function <- function(x, ..., args) {
  if (missing(args)) stop("'args' must be specified.")
  psis <- psislw(..., llfun = x, llargs = args, COMPUTE_LOOS = TRUE)
  out <- pointwise_loo(psis = psis, llfun = x, llargs = args)
  structure(out, log_lik_dim = with(args, c(S,N)), class = "loo")
}

