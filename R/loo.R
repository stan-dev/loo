#' Leave-one-out cross-validation (LOO)
#'
#' @description Efficient approximate leave-one-out cross-validation for
#'   Bayesian models using Pareto smoothed importance sampling (PSIS). See
#'   Vehtari, Gelman, and Gabry (2016, 2017) and \link{loo-package} for
#'   background.
#'
#'   The \code{loo} function is an S3 generic and methods are provided for
#'   computing LOO from 3-D pointwise log-likelihood arrays, pointwise
#'   log-likelihood matrices, and log-likelihood functions. The array and matrix
#'   methods are usually most convenient, but for models fit to very large
#'   datasets the \code{loo.function} method is more memory efficient and may be
#'   preferable.
#'
#'
#' @export loo loo.array loo.matrix loo.function
#' @param x A log-likelihood array, matrix, or function. See the \strong{Methods
#'   (by class)} section below for a detailed description of how to specify the
#'   inputs for each method.
#' @param ... Arguments passed on to the various methods.
#' @inheritParams psis
#'
#' @return A named list with class \code{c("psis_loo", "loo")} and components:
#' \describe{
#'  \item{\code{estimates}}{
#'   A matrix with two columns (\code{Estimate}, \code{SE}) and three rows
#'   (\code{elpd_loo}, \code{p_loo}, \code{looic}). This contains point
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
#' @seealso
#' \itemize{
#'  \item \code{\link{psis}} for the underlying Pareto Smoothed Importance
#'  Sampling (PSIS) procedure used in the LOO-CV approximation.
#'  \item \link{pareto-k-diagnostic} for convenience functions for looking at
#'  diagnostics.
#'  \item \code{\link{compare_models}} for model comparison.
#'  \item \code{\link{print.loo}} for info on the \code{print} method.
#' }
#' @template loo-and-psis-references
#'
#' @examples
#'
#' ### Array and matrix methods (using example objects included with loo package)
#' # Array method
#' LLarr <- example_loglik_array()
#' loo(LLarr)
#'
#' # Matrix method
#' LLmat <- example_loglik_matrix()
#' chain <- rep(1:ncol(LLarr), each = nrow(LLarr))
#' loo(LLmat, chain_id = chain)
#'
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
#' ### Using log-likelihood function instead of array or matrix
#' set.seed(024)
#'
#' # Simulate data and draw from posterior
#' N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
#' p <- rbeta(1, a0, b0)
#' y <- rbinom(N, size = K, prob = p)
#' a <- a0 + sum(y); b <- b0 + N * K - sum(y)
#' fake_posterior <- as.matrix(rbeta(S, a, b))
#' dim(fake_posterior) # S x 1
#' fake_data <- data.frame(y,K)
#' dim(fake_data) # N x 2
#'
#' llfun <- function(data_i, draws) {
#'   # each time called internally within loo the arguments will be equal to:
#'   # data_i: ith row of fake_data (fake_data[i,, drop=FALSE])
#'   # draws: entire fake_posterior matrix
#'   dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
#' }
#'
#' # Function method
#' loo_with_fn <- loo(
#'   x = llfun,
#'   chain_id = rep(1, nrow(fake_posterior)), # pretend all draws came from 1 chain for this example
#'   draws = fake_posterior,
#'   data = fake_data
#' )
#'
#' # Check that we get same answer if using log-likelihood matrix
#' mat <- sapply(1:N, function(i) {
#'   llfun(data_i = fake_data[i,, drop=FALSE], draws = fake_posterior)
#' })
#' loo_with_mat <- loo(mat, chain_id = rep(1, 100))
#' all.equal(loo_with_mat, loo_with_fn) # should be TRUE!
#'
loo <- function(x, ...) {
  UseMethod("loo")
}

#' @export
#' @templateVar fn loo
#' @template array
#'
loo.array <-
  function(x,
           ...,
           cores = getOption("loo.cores", 1)) {
    psis_out <- psis.array(-x, cores = cores)
    ll <- llarray_to_matrix(x)
    pointwise <- pointwise_loo_calcs(ll, psis_out)
    psis_loo_object(
      pointwise = pointwise,
      diagnostics = psis_out[c("pareto_k", "n_eff")],
      log_lik_dim = attr(psis_out, "dims")
    )
  }

#' @export
#' @templateVar fn loo
#' @template matrix
#' @template chain_id
#'
loo.matrix <-
  function(x,
           ...,
           chain_id,
           cores = getOption("loo.cores", 1)) {
    psis_out <- psis.matrix(-x, chain_id, cores = cores)
    pointwise <- pointwise_loo_calcs(x, psis_out)
    psis_loo_object(
      pointwise = pointwise,
      diagnostics = psis_out[c("pareto_k", "n_eff")],
      log_lik_dim = attr(psis_out, "dims")
    )
  }

#' @export
#' @templateVar fn loo
#' @template function
#' @param draws,data For the function method only. See the \strong{Methods (by
#'   class)} section below for details on these arguments.
#'
loo.function <-
  function(x,
           ...,
           chain_id,
           draws = NULL,
           data = NULL,
           cores = getOption("loo.cores", 1)) {

    stopifnot(is.data.frame(data) || is.matrix(data),
              !is.null(dim(draws)))
    S <- dim(draws)[1]
    N <- dim(data)[1]

    .llfun <- match.fun(x)
    arg_names <- names(formals(.llfun))
    if (length(arg_names) != 2 ||
        !all(arg_names %in% c("data_i", "draws"))) {
      stop("Log-likelihood function should have two arguments: ",
           "'data_i' and 'draws'.", call. = FALSE)
    }


    if (cores == 1) {
      psis_list <-
        lapply(
          X = seq_len(N),
          FUN = .loo_i,
          llfun = .llfun,
          data = data,
          draws = draws,
          chain_id = chain_id
        )
    } else {
      if (.Platform$OS.type != "windows") {
        psis_list <-
          parallel::mclapply(
            mc.cores = cores,
            X = seq_len(N),
            FUN = .loo_i,
            llfun = .llfun,
            data = data,
            draws = draws,
            chain_id = chain_id
          )
      } else {
        cl <- parallel::makePSOCKcluster(cores)
        on.exit(parallel::stopCluster(cl))
        psis_list <-
          parallel::parLapply(
            cl = cl,
            X = seq_len(N),
            fun = .loo_i,
            llfun = .llfun,
            data = data,
            draws = draws,
            chain_id = chain_id
          )
      }
    }

    pointwise <- do.call(rbind, lapply(psis_list, "[[", "estimates"))
    diagnostics <- list(
      pareto_k = psis_apply(psis_list, "pareto_k"),
      n_eff = psis_apply(psis_list, "n_eff")
    )

    psis_loo_object(
      pointwise = pointwise,
      diagnostics = diagnostics,
      log_lik_dim = c(S = S, N = N)
    )
}



# internal ----------------------------------------------------------------

#' Compute pointwise elpd_loo, p_loo, looic from log lik matrix and
#' psis log weights
#'
#' @noRd
#' @param ll Log-likelihood matrix.
#' @param psis_object The object returned by psis.
#' @return Named list with pointwise elpd_loo, p_loo, and looic.
#'
pointwise_loo_calcs <- function(ll, psis_object) {
  lpd <- logColMeansExp(ll)
  ll_plus_lw <- ll + weights(psis_object, normalize = TRUE, log = TRUE)
  elpd_loo <- colLogSumExps(ll_plus_lw)
  looic <- -2 * elpd_loo
  p_loo <- lpd - elpd_loo
  cbind(elpd_loo, p_loo, looic)
}

#' Structure the object returned by the loo methods
#'
#' @noRd
#' @param pointwise Matrix containing columns elpd_loo, p_loo, looic
#' @param diagnostics Named list containing vector 'pareto_k' and vector 'n_eff'
#' @param log_lik_dim Log likelihood matrix dimensions (attribute of psis object)
#'
psis_loo_object <- function(pointwise, diagnostics, log_lik_dim) {
  stopifnot(is.matrix(pointwise), is.list(diagnostics))
  estimates <- table_of_estimates(pointwise)
  structure(
    nlist(estimates, pointwise, diagnostics),
    log_lik_dim = log_lik_dim,
    class = c("psis_loo", "loo")
  )
}


#' @description The \code{loo_i} function enables testing log-likelihood
#'   functions for use with the \code{loo.function} method.
#'
#' @rdname loo
#' @export
#'
#' @param i For \code{loo_i}, an integer in 1:N.
#' @param llfun For \code{loo_i}, the same as \code{x} for the
#'   \code{loo.function} method. A log-likelihood function as described in the
#'   \strong{Methods (by class)} section.
#'
#' @return \code{loo_i} returns a named list of with results for a
#'   single observation. The list has the following components:
#'   \itemize{
#'     \item \code{estimates}: 1 x 3 matrix with columns \code{elpd_loo},
#'     \code{p_loo}, and \code{looic}.
#'     \item \code{pareto_k}: scalar \link[=pareto-k-diagnostic]{pareto k diagnostic}.
#'     \item \code{n_eff}: scalar \link[=psis]{PSIS} effective sample size estimate.
#'   }
#'
loo_i <-
  function(i,
           llfun,
           data,
           draws,
           chain_id = NULL) {
    stopifnot(
      i == as.integer(i),
      is.data.frame(data) || is.matrix(data),
      !is.null(dim(draws))
    )
    i <- as.integer(i)
    S <- dim(draws)[1]
    N <- dim(data)[1]
    stopifnot(i %in% seq_len(N))

    .loo_i(i, llfun = match.fun(llfun), data, draws, chain_id)
  }


# Function that is passed to the FUN argument of lapply, mclapply, or parLapply
# for the loo.function method. The arguments and return value are the same as
# the ones documented above for loo_i.
.loo_i <-
  function(i,
           llfun,
           data,
           draws,
           chain_id = NULL) {

    d_i <- data[i, , drop = FALSE]
    ll_i <- llfun(data_i = d_i, draws = draws)
    psis_out <- psis(x = -ll_i, chain_id = chain_id, cores = 1)
    list(
      estimates = pointwise_loo_calcs(ll_i, psis_out),
      pareto_k = psis_out$pareto_k,
      n_eff = psis_out$n_eff
    )
  }
