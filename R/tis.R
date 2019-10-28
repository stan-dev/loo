#' @rdname psis
#' @export
tis <- function(log_ratios, ...) UseMethod("tis")

#' @export
#' @templateVar fn psis
#' @template array
#'
tis.array <-
  function(log_ratios, ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
  importance_sampling.array(log_ratios = log_ratios, ...,
                            r_eff = r_eff,
                            cores = cores,
                            is_method = "tis")
  }

#' @export
#' @templateVar fn psis
#' @template matrix
#'
tis.matrix <-
  function(log_ratios,
           ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
    importance_sampling.matrix(log_ratios,
                               ...,
                               r_eff = r_eff,
                               cores = cores,
                               is_method = "tis")
  }

#' @export
#' @templateVar fn psis
#' @template vector
#'
tis.default <-
  function(log_ratios, ..., r_eff = NULL) {
    importance_sampling.default(log_ratios = log_ratios, ...,
                                r_eff = r_eff, is_method = "tis")
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
#' @param tail_len_i Not used. Included to conform to PSIS API.
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
do_tis_i <- function(log_ratios_i, tail_len_i) {
  S <- length(log_ratios_i)
  log_Z <- logSumExp(log_ratios_i) - log(S) # Normalization term, c-hat in Ionides (2008) appendix
  log_cutpoint <- log_Z + 0.5 * log(S)
  lw_i <- pmin(log_ratios_i, log_cutpoint)
  lw_i <- lw_i - max(lw_i)
  list(log_weights = lw_i, pareto_k = 0)
}

