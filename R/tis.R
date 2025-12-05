#' Truncated importance sampling (TIS)
#'
#' Implementation of truncated (self-normalized) importance sampling (TIS),
#' truncated at S^(1/2) as recommended by Ionides (2008).
#'
#' @param log_ratios An array, matrix, or vector of importance ratios on the log
#'   scale (for Importance sampling LOO, these are *negative* log-likelihood
#'   values). See the **Methods (by class)** section below for a detailed
#'   description of how to specify the inputs for each method.
#' @template cores
#' @param ... Arguments passed on to the various methods.
#' @param r_eff Vector of relative effective sample size estimates containing
#'   one element per observation. The values provided should be the relative
#'   effective sample sizes of `1/exp(log_ratios)` (i.e., `1/ratios`).
#'   This is related to the relative efficiency of estimating the normalizing
#'   term in self-normalizing importance sampling. If `r_eff` is not
#'   provided then the reported (T)IS effective sample sizes and Monte Carlo
#'   error estimates can be over-optimistic. If the posterior draws are (near)
#'   independent then `r_eff=1` can be used. `r_eff` has to be a scalar (same
#'   value is used for all observations) or a vector with length equal to the
#'   number of observations. The default value is 1. See the [relative_eff()]
#'   helper function for computing `r_eff`.
#'
#' @return The `tis()` methods return an object of class `"tis"`,
#'   which is a named list with the following components:
#'
#' \describe{
#'   \item{`log_weights`}{
#'     Vector or matrix of smoothed (and truncated) but *unnormalized* log
#'     weights. To get normalized weights use the
#'     [`weights()`][weights.importance_sampling] method provided for objects of
#'     class `tis`.
#'   }
#'  \item{`diagnostics`}{
#'    A named list containing one vector:
#'    * `pareto_k`: Not used in `tis`, all set to 0.
#'    * `n_eff`: Effective sample size estimates.
#'  }
#' }
#'
#' Objects of class `"tis"` also have the following [attributes][attributes()]:
#' \describe{
#'   \item{`norm_const_log`}{
#'     Vector of precomputed values of `colLogSumExps(log_weights)` that are
#'     used internally by the [weights()]method to normalize the log weights.
#'   }
#'   \item{`r_eff`}{
#'     If specified, the user's `r_eff` argument.
#'   }
#'   \item{`tail_len`}{
#'     Not used for `tis`.
#'   }
#'   \item{`dims`}{
#'     Integer vector of length 2 containing `S` (posterior sample size)
#'     and `N` (number of observations).
#'   }
#'   \item{`method`}{
#'     Method used for importance sampling, here `tis`.
#'   }
#' }
#'
#' @seealso
#' * [psis()] for approximate LOO-CV using PSIS.
#' * [loo()] for approximate LOO-CV.
#' * [pareto-k-diagnostic] for PSIS diagnostics.
#'
#' @references
#' Ionides, Edward L. (2008). Truncated importance sampling.
#' *Journal of Computational and Graphical Statistics* 17(2): 295--311.
#'
#' @examples
#' log_ratios <- -1 * example_loglik_array()
#' r_eff <- relative_eff(exp(-log_ratios))
#' tis_result <- tis(log_ratios, r_eff = r_eff)
#' str(tis_result)
#'
#' # extract smoothed weights
#' lw <- weights(tis_result) # default args are log=TRUE, normalize=TRUE
#' ulw <- weights(tis_result, normalize=FALSE) # unnormalized log-weights
#'
#' w <- weights(tis_result, log=FALSE) # normalized weights (not log-weights)
#' uw <- weights(tis_result, log=FALSE, normalize = FALSE) # unnormalized weights
#'
#' @export
tis <- function(log_ratios, ...) UseMethod("tis")

#' @export
#' @templateVar fn tis
#' @template array
#'
tis.array <-
  function(log_ratios, ...,
           r_eff = 1,
           cores = getOption("mc.cores", 1)) {
  importance_sampling.array(log_ratios = log_ratios, ...,
                            r_eff = r_eff,
                            cores = cores,
                            method = "tis")
  }

#' @export
#' @templateVar fn tis
#' @template matrix
#'
tis.matrix <-
  function(log_ratios,
           ...,
           r_eff = 1,
           cores = getOption("mc.cores", 1)) {
    importance_sampling.matrix(log_ratios,
                               ...,
                               r_eff = r_eff,
                               cores = cores,
                               method = "tis")
  }

#' @export
#' @templateVar fn tis
#' @template vector
#'
tis.default <-
  function(log_ratios, ..., r_eff = 1) {
    importance_sampling.default(log_ratios = log_ratios, ...,
                                r_eff = r_eff, method = "tis")
  }


#' @rdname psis
#' @export
is.tis <- function(x) {
  inherits(x, "tis") && is.list(x)
}


# internal ----------------------------------------------------------------

#' Truncated Importance Sampling on a single vector
#'
#' @noRd
#' @param log_ratios_i A vector of log importance ratios (for `loo()`, negative
#'   log likelihoods).
#' @param ... Not used. Included to conform to PSIS API.
#'
#' @details Implementation of Truncated importance sampling (TIS), a method for
#' stabilizing importance ratios. The version of TIS implemented here
#' corresponds to the algorithm presented in Ionides (2008) with truncation at
#' sqrt(S).
#'
#' @return A named list containing:
#' * `lw`: vector of unnormalized log weights
#' * `pareto_k`: scalar Pareto k estimate. For 'tis', this defaults  to 0.
#'
do_tis_i <- function(log_ratios_i, ...) {
  S <- length(log_ratios_i)
  log_Z <- logSumExp(log_ratios_i) - log(S) # Normalization term, c-hat in Ionides (2008) appendix
  log_cutpoint <- log_Z + 0.5 * log(S)
  lw_i <- pmin(log_ratios_i, log_cutpoint)
  list(log_weights = lw_i, pareto_k = 0)
}

