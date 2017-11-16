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
#' @param r_eff Vector of relative effective sample size estimates for the
#'   likelihood (\code{exp(log_lik)}) of each observation. This is related to
#'   the relative efficiency of estimating the normalizing term in
#'   self-normalizing importance sampling. The default is \code{NULL}, in which
#'   case Monte Carlo error estimates are not computed. See the
#'   \code{\link{relative_eff}} helper function for computing \code{r_eff}.
#' @template cores
#'
#' @return A named list with class \code{c("psis_loo", "loo")} and components:
#' \describe{
#'  \item{\code{estimates}}{
#'   A matrix with two columns (\code{Estimate}, \code{SE}) and four rows
#'   (\code{elpd_loo}, \code{mcse_elpd_loo}, \code{p_loo}, \code{looic}). This
#'   contains point estimates and standard errors of the expected log pointwise
#'   predictive density (\code{elpd_loo}), the Monte Carlo standard error of
#'   \code{elpd_loo} (\code{mcse_elpd_loo}), the effective number of parameters
#'   (\code{p_loo}) and the LOO information criterion \code{looic} (which is
#'   just \code{-2 * elpd_loo}, i.e., converted to deviance scale).
#'  }
#'  \item{\code{pointwise}}{
#'   A matrix with four columns (and number of rows equal to the number of
#'   observations) containing the pointwise contributions of each of the above
#'   measures (\code{elpd_loo}, \code{mcse_elpd_loo}, \code{p_loo},
#'   \code{looic}).
#'  }
#'  \item{\code{diagnostics}}{
#'  A named list containing two vectors:
#'   \itemize{
#'    \item \code{pareto_k}: Estimates of the shape parameter \eqn{k} of the
#'    generalized Pareto fit to the importance ratios for each leave-one-out
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
#' rel_n_eff <- relative_eff(exp(LLarr))
#' loo(LLarr), r_eff = rel_n_eff, cores = 2)
#'
#' # Matrix method
#' LLmat <- example_loglik_matrix()
#' rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = 500))
#' loo(LLmat, r_eff = rel_n_eff, cores = 2)
#'
#'
#' \dontrun{
#' ### Usage with stanfit objects
#' # see ?extract_log_lik
#' log_lik1 <- extract_log_lik(stanfit1, merge_chains = FALSE)
#' rel_n_eff <- relative_eff(exp(log_lik1))
#' loo(log_lik1, r_eff = rel_n_eff, cores = 2)
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
#'   draws = fake_posterior,
#'   data = fake_data
#' )
#'
#' # Check that we get same answer if using log-likelihood matrix
#' mat <- sapply(1:N, function(i) {
#'   llfun(data_i = fake_data[i,, drop=FALSE], draws = fake_posterior)
#' })
#' loo_with_mat <- loo(mat)
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
           r_eff = NULL,
           cores = getOption("loo.cores", 1),
           save_psis = FALSE) {
    throw_r_eff_warning(r_eff)
    psis_out <- psis.array(log_ratios = -x, r_eff = r_eff, cores = cores)
    ll <- llarray_to_matrix(x)
    pointwise <- pointwise_loo_calcs(ll, psis_out)
    psis_loo_object(
      pointwise = pointwise,
      diagnostics = psis_out$diagnostics,
      dims = dim(psis_out),
      psis_object = if (save_psis) psis_out else NULL
    )
  }

#' @export
#' @templateVar fn loo
#' @template matrix
#'
loo.matrix <-
  function(x,
           ...,
           r_eff = NULL,
           cores = getOption("loo.cores", 1),
           save_psis = FALSE) {
    throw_r_eff_warning(r_eff)
    psis_out <- psis.matrix(log_ratios = -x, r_eff = r_eff, cores = cores)
    pointwise <- pointwise_loo_calcs(x, psis_out)
    psis_loo_object(
      pointwise = pointwise,
      diagnostics = psis_out$diagnostics,
      dims = dim(psis_out),
      psis_object = if (save_psis) psis_out else NULL
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
           data = NULL,
           draws = NULL,
           r_eff = NULL,
           cores = getOption("loo.cores", 1)) {

    stopifnot(is.data.frame(data) || is.matrix(data),
              !is.null(dim(draws)))
    throw_r_eff_warning(r_eff)
    S <- dim(draws)[1]
    N <- dim(data)[1]

    .llfun <- validate_llfun(x)

    if (cores == 1) {
      psis_list <-
        lapply(
          X = seq_len(N),
          FUN = .loo_i,
          llfun = .llfun,
          data = data,
          draws = draws,
          r_eff = r_eff
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
            r_eff = r_eff
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
            r_eff = r_eff
          )
      }
    }

    pointwise <- lapply(psis_list, "[[", "estimates")
    diagnostics <- lapply(psis_list, "[[", "diagnostics")

    psis_loo_object(
      pointwise = do.call(rbind, pointwise),
      diagnostics = list(
        pareto_k = psis_apply(diagnostics, "pareto_k"),
        n_eff = psis_apply(diagnostics, "n_eff")
      ),
      dims = c(S, N)
    )
  }


#' @description The \code{loo_i} function enables testing log-likelihood
#'   functions for use with the \code{loo.function} method.
#'
#' @rdname loo
#' @export
#'
#' @param i For \code{loo_i}, an integer in \code{1:N}.
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
           r_eff = NULL) {
    stopifnot(
      i == as.integer(i),
      is.data.frame(data) || is.matrix(data),
      !is.null(dim(draws))
    )
    i <- as.integer(i)
    S <- dim(draws)[1]
    N <- dim(data)[1]
    stopifnot(i %in% seq_len(N))

    .loo_i(
      i = i,
      llfun = match.fun(llfun),
      data = data,
      draws = draws,
      r_eff = r_eff[i]
    )
  }

#' @export
dim.psis_loo <- function(x) {
  attr(x, "dims")
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
  mcse_elpd_loo <- mcse_elpd(ll, elpd_loo, psis_object)
  cbind(elpd_loo, mcse_elpd_loo, p_loo, looic)
}

#' Structure the object returned by the loo methods
#'
#' @noRd
#' @param pointwise Matrix containing columns elpd_loo, p_loo, looic
#' @param diagnostics Named list containing vector 'pareto_k' and vector 'n_eff'
#' @param dims Log likelihood matrix dimensions (attribute of psis object)
#' @param psis_object PSIS object.
#'
psis_loo_object <- function(pointwise, diagnostics, dims, psis_object = NULL) {
  stopifnot(is.matrix(pointwise), is.list(diagnostics))
  estimates <- table_of_estimates(pointwise)
  structure(
    nlist(estimates, pointwise, diagnostics, psis_object),
    dims = dims,
    class = c("psis_loo", "loo")
  )
}


# Function that is passed to the FUN argument of lapply, mclapply, or parLapply
# for the loo.function method. The arguments and return value are the same as
# the ones documented above for loo_i.
.loo_i <-
  function(i,
           llfun,
           data,
           draws,
           r_eff = NULL) {

    if (!is.null(r_eff)) {
      r_eff <- r_eff[i]
    }
    d_i <- data[i, , drop = FALSE]
    ll_i <- llfun(data_i = d_i, draws = draws)
    psis_out <- psis(log_ratios = -ll_i, r_eff = r_eff, cores = 1)
    list(
      estimates = pointwise_loo_calcs(ll_i, psis_out),
      diagnostics = psis_out$diagnostics
    )
  }


#' Compute Monte Carlo standard error for ELPD
#'
#' @noRd
#' @param ll Log-likelihood matrix.
#' @param E_elpd elpd_loo column of pointwise matrix.
#' @param psis_object Object returned by psis.
#' @param n_samples Number of draws to take from N(E[epd_i], SD[epd_i]).
#' @return Vector of standard error estimates
#'
mcse_elpd <- function(ll, E_elpd, psis_object, n_samples = 1000) {
  E_epd <- exp(E_elpd)
  w <- weights(psis_object, log = FALSE)
  r_eff <- attr(psis_object, "r_eff")
  var_elpd <- sapply(seq_len(ncol(w)), function(i) {
    var_epd_i <- sum(w[, i]^2 * (exp(ll[, i]) - E_epd[i])^2)
    sd_epd_i <- sqrt(var_epd_i)
    z <- rnorm(n_samples, mean = E_epd[i], sd = sd_epd_i)
    var(log(z))
  })
  sqrt(var_elpd / r_eff)
}


#' Warning message if r_eff not specified
#'
#' @noRd
#' @param r_eff User's r_eff argument or NULL.
#' @return Nothing, just throws a warning if r_eff is NULL.
#'
throw_r_eff_warning <- function(r_eff = NULL) {
  if (is.null(r_eff)) {
    warning(
      "Relative effective sample sizes ('r_eff' argument) not specified.\n",
      "For models fit with MCMC, the reported PSIS effective sample sizes and \n",
      "MCSE estimates will be over-optimistic.",
      call. = FALSE
    )
  }
}
