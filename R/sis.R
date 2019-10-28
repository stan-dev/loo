#' @rdname psis
#' @export
sis <- function(log_ratios, ...) UseMethod("sis")

#' @export
#' @templateVar fn psis
#' @template array
#'
sis.array <-
  function(log_ratios, ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
  importance_sampling.array(log_ratios = log_ratios, ...,
                            r_eff = r_eff,
                            cores = cores,
                            is_method = "sis")
  }

#' @export
#' @templateVar fn psis
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
                               is_method = "sis")
  }

#' @export
#' @templateVar fn psis
#' @template vector
#'
sis.default <-
  function(log_ratios, ..., r_eff = NULL) {
    importance_sampling.default(log_ratios = log_ratios, ...,
                                r_eff = r_eff, is_method = "sis")
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
#' @param tail_len_i Not used. Included to conform to PSIS API.
#'
#' @details Implementation standard importance sampling.
#' @return A named list containing:
#' * `lw`: vector of unnormalized log weights
#' * `pareto_k`: scalar Pareto k estimate. For IS, this defaults to 0.
do_sis_i <- function(log_ratios_i, tail_len_i) {
  S <- length(log_ratios_i)
  lw_i <- log_ratios_i - max(log_ratios_i)
  list(log_weights = lw_i, pareto_k = 0)
}
