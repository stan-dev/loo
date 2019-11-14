#' Standard importance sampling (SIS)
#'
#' Implementation of standard importance sampling (SIS).
#'
#' @param log_ratios An array, matrix, or vector of importance ratios on the log
#'   scale (for Importance sampling LOO, these are *negative* log-likelihood values). See the
#'   **Methods (by class)** section below for a detailed description of how
#'   to specify the inputs for each method.
#' @template cores
#' @param ... Arguments passed on to the various methods.
#' @param r_eff Vector of relative effective sample size estimates containing
#'   one element per observation. The values provided should be the relative
#'   effective sample sizes of `1/exp(log_ratios)` (i.e., `1/ratios`).
#'   This is related to the relative efficiency of estimating the normalizing
#'   term in self-normalizing importance sampling. See the [relative_eff()]
#'   helper function for computing `r_eff`. If using `psis` with
#'   draws of the `log_ratios` not obtained from MCMC then the warning
#'   message thrown when not specifying `r_eff` can be disabled by
#'   setting `r_eff` to `NA`.
#'
#' @return The `sis()` methods return an object of class `"sis"`,
#'   which is a named list with the following components:
#'
#' \describe{
#'   \item{`log_weights`}{
#'     Vector or matrix of smoothed (and truncated) but *unnormalized* log
#'     weights, *minus the largest log ratio* for numerical reasons.
#'     To get normalized weights use the `weights` method provided
#'     for objects of class `sis`.
#'   }
#'  \item{`diagnostics`}{
#'    A named list containing one vector:
#'    * `pareto_k`: Not used in `sis`, all set to 0.
#'    * `n_eff`: effective sample size estimates.
#'  }
#' }
#'
#' Objects of class `"sis"` also have the following [attributes][attributes()]:
#' \describe{
#'   \item{`norm_const_log`}{
#'     Vector of precomputed values of `colLogSumExps(log_weights)` that are
#'     used internally by the `weights` method to normalize the log weights.
#'   }
#'   \item{`r_eff`}{
#'     If specified, the user's `r_eff` argument.
#'   }
#'   \item{`tail_len`}{
#'     Not used for `sis`.
#'   }
#'   \item{`dims`}{
#'     Integer vector of length 2 containing `S` (posterior sample size)
#'     and `N` (number of observations).
#'   }
#'   \item{`method`}{
#'     Method used for importance sampling, here `sis`.
#'   }
#' }
#'
#' @seealso
#' \itemize{
#' \item [psis()] for approximate LOO-CV using PSIS.
#' \item [loo()] for approximate LOO-CV.
#' \item [pareto-k-diagnostic] for PSIS diagnostics.
#' }
#'
#' @template loo-and-psis-references
#'
#' @examples
#' log_ratios <- -1 * example_loglik_array()
#' r_eff <- relative_eff(exp(-log_ratios))
#' sis_result <- sis(log_ratios, r_eff = r_eff)
#' str(sis_result)
#'
#' # extract smoothed weights
#' lw <- weights(sis_result) # default args are log=TRUE, normalize=TRUE
#' ulw <- weights(sis_result, normalize=FALSE) # unnormalized log-weights
#'
#' w <- weights(sis_result, log=FALSE) # normalized weights (not log-weights)
#' uw <- weights(sis_result, log=FALSE, normalize = FALSE) # unnormalized weights
#'
#' @export
sis <- function(log_ratios, ...) UseMethod("sis")

#' @export
#' @templateVar fn sis
#' @template array
#'
sis.array <-
  function(log_ratios, ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
  importance_sampling.array(log_ratios = log_ratios, ...,
                            r_eff = r_eff,
                            cores = cores,
                            method = "sis")
  }

#' @export
#' @templateVar fn sis
#' @template matrix
#'
sis.matrix <-
  function(log_ratios,
           ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
    importance_sampling.matrix(log_ratios,
                               ...,
                               r_eff = r_eff,
                               cores = cores,
                               method = "sis")
  }

#' @export
#' @templateVar fn sis
#' @template vector
#'
sis.default <-
  function(log_ratios, ..., r_eff = NULL) {
    importance_sampling.default(log_ratios = log_ratios, ...,
                                r_eff = r_eff,
                                method = "sis")
  }

#' @rdname psis
#' @export
is.sis <- function(x) {
  inherits(x, "sis") && is.list(x)
}


# internal ----------------------------------------------------------------


#' Standard IS on a single vector
#'
#' @noRd
#' @param log_ratios_i A vector of log importance ratios (for `loo()`, negative
#'   log likelihoods).
#' @param ... Not used. Included to conform to PSIS API.
#'
#' @details Implementation standard importance sampling.
#' @return A named list containing:
#' * `lw`: vector of unnormalized log weights
#' * `pareto_k`: scalar Pareto k estimate. For IS, this defaults to 0.
do_sis_i <- function(log_ratios_i, ...) {
  S <- length(log_ratios_i)
  lw_i <- log_ratios_i - max(log_ratios_i)
  list(log_weights = lw_i, pareto_k = 0)
}
