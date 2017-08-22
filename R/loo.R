#' Leave-one-out cross-validation (LOO)
#'
#' Efficient approximate leave-one-out cross-validation for Bayesian models. See
#' \link{loo-package} and Vehtari, Gelman, and Gabry (2016, 2017) for
#' background.
#'
#' @export loo loo.array loo.matrix loo.function
#' @inheritParams psis
#' @param x A log-likelihood array, matrix, or function. See the \strong{Methods
#'   (by class)} section below for a detailed description of how to specify the
#'   inputs for each method.
#' @param ... Additional arguments passed on to \code{\link{psis}}. Possible
#'   arguments and their defaults are:
#' \describe{
#'  \item{\code{wtrunc = 3/4}}{
#'   For truncating very large weights to \eqn{S}^\code{wtrunc} (set to zero
#'   for no truncation). We recommend the default value unless there are
#'   problems.
#'  }
#'}
#'
#' @return A named list with class \code{c("psis_loo", "loo")} and components:
#' \describe{
#'  \item{\code{estimates}}{
#'   A matrix with two columns (\code{"Estimate"}, \code{"SE"}) and three rows
#'   (\code{"elpd_loo"}, \code{"p_loo"}, \code{"looic"}). This contains point
#'   estimates and standard errors of the expected log pointwise predictive
#'   density (\code{elpd_loo}), the effective number of parameters
#'   (\code{p_loo}) and the LOO information criterion \code{looic} (which is
#'   just \code{-2 * elpd_loo}, i.e., converted to deviance scale).
#'  }
#'  \item{\code{pointwise}}{
#'   A matrix with three columns (and number of rows equal to the number of
#'   observations) containing the pointwise contributions of each of the above
#'   measures (\code{elpd_loo}, \code{p_loo}, \code{looic}).
#'  }
#'  \item{\code{diagnostics}}{
#'  A named list containing two vectors:
#'   \itemize{
#'    \item \code{pareto_k}: Estimates of the shape parameter \eqn{k} of the
#'    generaelized Pareto fit to the importance ratios for each leave-one-out
#'    distribution. See the \code{\link{pareto-k-diagnostic}} page for details.
#'    \item \code{n_eff}: PSIS effective sample size estimates.
#'   }
#'  }
#' }
#'
#' @note For models fit to very large datasets we recommend the
#'   \code{loo.function} method, which is much more memory efficient than the
#'   array and matrix methods. However, the array and matrix methods are
#'   typically more convenient, so it is usually worth trying them and then
#'   switching to \code{loo.function} if memory is an issue.
#'
#' @seealso
#' \itemize{
#'  \item \code{\link{psis}} for the underlying Pareto Smoothed Importance
#'  Sampling (PSIS) procedure used in the LOO-CV approximation.
#'  \item \link{pareto-k-diagnostic} for convenience functions for looking at
#'  diagnostics.
#'  \item \code{\link{compare}} for model comparison.
#'  \item \code{\link{print.loo}} for info on the \code{print} method.
#' }
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
#' llfun <- function(data, draws) {
#'   dbinom(data$y, size = data$K, prob = draws, log = TRUE)
#' }
#' loo_with_fn <- loo(llfun, chain_id = rep(1, 100), args = nlist(data, draws, N, S), cores = 1)
#'
#' # Check that we get same answer if using log-likelihood matrix
#' log_lik_mat <- sapply(1:N, function(i) llfun(data[i,, drop=FALSE], draws))
#' loo_with_mat <- loo(log_lik_mat, chain_id = rep(1, 100), cores = 1)
#' all.equal(loo_with_mat, loo_with_fn)
#' }
#'
loo <- function(x, ...) {
  UseMethod("loo")
}

#' @export
#' @templateVar fn loo
#' @template array
#'
loo.array <- function(x, ..., cores = getOption("loo.cores", 1)) {
  psis_out <- psis.array(x, ..., cores = cores)
  ll <- llarray_to_matrix(x)
  pointwise <- pointwise_loo_calcs(ll, psis_out)
  psis_loo_object(
    pointwise = pointwise,
    diagnostics = psis_out[c("pareto_k", "n_eff")],
    log_lik_dim = attr(psis_out, "log_lik_dim")
  )
}

#' @export
#' @templateVar fn loo
#' @template matrix
#' @template chain_id
#'
loo.matrix <- function(x, ..., chain_id, cores = getOption("loo.cores", 1)) {
  psis_out <- psis.matrix(x, chain_id, ...)
  pointwise <- pointwise_loo_calcs(x, psis_out)
  psis_loo_object(
    pointwise = pointwise,
    diagnostics = psis_out[c("pareto_k", "n_eff")],
    log_lik_dim = attr(psis_out, "log_lik_dim")
  )
}

#' @export
#' @templateVar fn loo
#' @template function
#'
loo.function <-
  function(x,
           ...,
           chain_id,
           draws = NULL,
           data = data.frame(),
           S = integer(),
           N = integer(),
           cores = getOption("loo.cores", 1)) {

    stopifnot(is.data.frame(data) || is.matrix(data))
    .LogLik <- match.fun(x)

    if (cores == 1) {
      psis_list <- lapply(
        X = seq_len(N),
        FUN = function(i) {
          d_i <- data[i,, drop=FALSE]
          ll_i <- as.matrix(.LogLik(d_i, draws))
          psis_out <- psis.matrix(as.matrix(ll_i), chain_id, cores = 1, ...)
          list(
            pointwise = pointwise_loo_calcs(ll_i, psis_out),
            diagnostics = cbind(psis_out$pareto_k, psis_out$n_eff)
          )
        }
      )
      # for (i in 1:N) {
      #   d_i <- data[i,, drop=FALSE]
      #   ll_i <- as.matrix(.LogLik(d_i, draws))
      #   psis_out <- psis.matrix(ll_i, chain_id, cores = 1, ...)
      #   pointwise_list[[i]] <- pointwise_loo_calcs(ll_i, psis_out)
      #   diagnostics_list[[i]] <- cbind(psis_out[["pareto_k"]], psis_out[["n_eff"]])
      # }
    } else {
      if (.Platform$OS.type != "windows") {
        psis_list <-
          parallel::mclapply(
            X = seq_len(N),
            FUN = function(i) {
              d_i <- data[i,, drop=FALSE]
              ll_i <- as.matrix(.LogLik(d_i, draws))
              psis_out <- psis.matrix(as.matrix(ll_i), chain_id, cores = 1, ...)
              list(
                pointwise = pointwise_loo_calcs(ll_i, psis_out),
                diagnostics = cbind(psis_out$pareto_k, psis_out$n_eff)
              )
            },
            mc.cores = cores
          )
      } else {
        cl <- parallel::makePSOCKcluster(cores)
        on.exit(parallel::stopCluster(cl))
        psis_list <-
          parallel::parLapply(
            cl = cl,
            X = seq_len(N),
            fun = function(i) {
              d_i <- data[i,, drop=FALSE]
              ll_i <- as.matrix(.LogLik(d_i, draws))
              psis_out <- psis.matrix(as.matrix(ll_i), chain_id, cores = 1, ...)
              list(
                pointwise = pointwise_loo_calcs(ll_i, psis_out),
                diagnostics = cbind(psis_out$pareto_k, psis_out$n_eff)
              )
            }
          )
      }
    }

    pointwise_list <- lapply(psis_list, "[[", "pointwise")
    diagnostics_list <- lapply(psis_list, "[[", "diagnostics")

    pointwise <- do.call(rbind, pointwise_list)
    diagnostics <- do.call(rbind, diagnostics_list)
    psis_loo_object(
      pointwise = pointwise,
      diagnostics = list(
        pareto_k = diagnostics[, 1],
        n_eff = diagnostics[, 2]
      ),
      log_lik_dim = c(S = S, N = N)
    )
}



# internal ----------------------------------------------------------------

# Compute pointwise elpd_loo, p_loo, looic from log lik matrix and
# psis log weights
#
# @param ll Log-likelihood matrix.
# @param psis_object The object returned by psis.
# @return Named list with pointwise elpd_loo, p_loo, and looic.
#
pointwise_loo_calcs <- function(ll, psis_object) {
  lpd <- logColMeansExp(ll)
  ll_plus_lw <- ll + weights(psis_object, normalize = TRUE, log = TRUE)
  elpd_loo <- colLogSumExps(ll_plus_lw)
  looic <- -2 * elpd_loo
  p_loo <- lpd - elpd_loo
  cbind(elpd_loo, p_loo, looic)
}

# structure the object returned by the loo methods
#
# @param pointwise Matrix containing columns elpd_loo, p_loo, looic
# @param diagnostics Named list containing vector 'pareto_k' and vector 'n_eff'
# @param log_lik_dim Log likelihood matrix dimensions (attribute of psis object)
#
psis_loo_object <- function(pointwise, diagnostics, log_lik_dim) {
  stopifnot(is.matrix(pointwise), is.list(diagnostics))
  estimates <- table_of_estimates(pointwise)
  structure(
    nlist(estimates, pointwise, diagnostics),
    log_lik_dim = log_lik_dim,
    class = c("psis_loo", "loo")
  )
}
