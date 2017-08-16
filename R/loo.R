#' Leave-one-out cross-validation (LOO)
#'
#' Efficient approximate leave-one-out cross-validation for Bayesian models. See
#' \link{loo-package} and Vehtari, Gelman, and Gabry (2016, 2017) for
#' background.
#'
#' @export loo loo.array loo.matrix loo.function
#' @param x A log-likelihood array, matrix, or function. See the \strong{Methods
#'   (by class)} section below for a detailed description.
#' @param args Only required if \code{x} is a function. A list containing
#'   the data required to specify the arguments to the function. See the
#'   \strong{Methods (by class)} section below for how \code{args} should be
#'   specified.
#' @param ... Optional arguments to pass to \code{\link{psis}}. Possible
#'   arguments and their defaults are:
#' \describe{
#'  \item{\code{wtrunc = 3/4}}{For truncating very large weights to
#'    \eqn{S}^\code{wtrunc} (set to zero for no truncation).
#'    We recommend the default value unless there are problems.}
#'  \item{\code{cores = getOption("loo.cores", 1)}}{
#'   The number of cores to use for parallelization. The default for
#'   an entire R session can be set with \code{options(loo.cores = NUMBER)}. As
#'   of \pkg{loo} version \code{2.0.0}, \strong{the default is 1 core}, but we
#'   recommend using as many (or close to as many) cores as possible.
#'  }
#'}
#'
#' @return A named list with class \code{c("psis_loo", "loo")} and components:
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
#' @note For models fit to very large datasets we recommend the
#'   \code{loo.function} method, which is much more memory efficient than the
#'   array and matrix methods. However, the array and matrix methods are
#'   typically more convenient, so it is usually worth trying them and then
#'   switching to \code{loo.function} if memory is an issue.
#'
#' @seealso
#' \code{\link{psis}} for the underlying Pareto Smoothed Importance Sampling
#' (PSIS) procedure used for approximating LOO.
#'
#' \code{\link{pareto-k-diagnostic}} for convenience functions for looking at
#' diagnostics.
#'
#' \code{\link{compare}} for model comparison.
#'
#' \code{\link{print.loo}} for the \code{print} method for \code{'loo'} objects.
#'
#' @template loo-and-psis-references
#'
#' @examples
#' \dontrun{
#' ### Usage with stanfit objects
#' # see ?extract_log_lik
#' log_lik1 <- extract_log_lik(stanfit1, merge_chains = FALSE)
#' loo1 <- loo(log_lik1)
#' print(loo1, digits = 3)
#'
#' log_lik2 <- extract_log_lik(stanfit2, merge_chains = FALSE)
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
#' @template array
#'
loo.array <- function(x, ...) {
  psis_out <- psis.array(x, ...)
  pointwise <- pointwise_loo_calcs(llarray_to_matrix(x), psis_out)
  nms <- c(names(pointwise), paste0("se_", names(pointwise)))
  psis_loo_object(
    stats = setNames(totals(pointwise), nms),
    pointwise = pointwise,
    diagnostics = psis_out[c("pareto_k", "n_eff")],
    lldim = attr(psis_out, "log_lik_dim")
  )
}

#' @export
#' @templateVar fn loo
#' @template matrix
#' @template chain_id
#'
loo.matrix <- function(x, chain_id, ...) {
  psis_out <- psis.matrix(x, chain_id, ...)
  pointwise <- pointwise_loo_calcs(x, psis_out)
  nms <- c(names(pointwise), paste0("se_", names(pointwise)))
  psis_loo_object(
    stats = setNames(totals(pointwise), nms),
    pointwise = pointwise,
    diagnostics = psis_out[c("pareto_k", "n_eff")],
    lldim = attr(psis_out, "log_lik_dim")
  )
}

#' @export
#' @templateVar fn loo
#' @template function
#'
loo.function <- function(x, args, ...) {
  psis <- psislw(..., llfun = x, llargs = args, COMPUTE_LOOS = TRUE)
  out <- pointwise_loo(psis = psis, llfun = x, llargs = args)
  structure(out, log_lik_dim = with(args, c(S,N)), class = c("psis_loo", "loo"))
}



# internal ----------------------------------------------------------------

# compute elpd_loo from log lik matrix and psis log weights
# @param ll log-likelihood matrix
# @param psis_object object returned by psis
pointwise_loo_calcs <- function(ll, psis_object) {
  lpd <- logColMeansExp(ll)
  ll_plus_lw <- ll + weights(psis_object, normalize = TRUE, log = TRUE)
  elpd_loo <- matrixStats::colLogSumExps(ll_plus_lw)
  looic <- -2 * elpd_loo
  p_loo <- lpd - elpd_loo
  nlist(elpd_loo, p_loo, looic)
}

# structure the object returned by the loo methods
# @param stats named list containing scalar values elpd_loo, p_loo, looic,
#   se_elpd_loo, se_p_loo, se_looic
# @param pointwise named list containing vectors elpd_loo, p_loo, looic
# @param diagnostics named list containing vector pareto_k and vector n_eff
#
psis_loo_object <- function(stats, pointwise, diagnostics, lldim) {
  out <- c(stats, diagnostics)
  out$pointwise <- do.call(cbind, pointwise)
  structure(
    out,
    log_lik_dim = lldim,
    class = c("psis_loo", "loo")
  )
}

