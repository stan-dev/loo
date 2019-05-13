#' Diagnostics for Laplace and ADVI approximations and Laplace-loo and ADVI-loo
#'
#' @param log_p The log-posterior (target) evaluated at S samples from the
#'   proposal distribution (q). A vector of length S.
#' @param log_q The log-density (proposal) evaluated at S samples from the
#'   proposal distribution (q). A vector of length S.
#' @param log_liks A log-likelihood matrix of size S * N, where N is the number
#'   of observations and S is the number of samples from q. See
#'   [loo.matrix()] for details. Default is `NULL`. Then only the
#'   posterior is evaluated using the k_hat diagnostic.
#' @inheritParams  loo
#'
#' @return
#' If log likelihoods are supplied, the function returns a `"loo"` object,
#' otherwise the function returns a `"psis"` object.
#'
#' @seealso [loo()] and [psis()]
#'
#' @template loo-and-psis-references
#'
#' @keywords internal
#'
psis_approximate_posterior <- function(log_p, log_q, log_liks = NULL, cores, save_psis, ...){
  if (!requireNamespace("checkmate", quietly=TRUE)) {
    stop("Please install the 'checkmate' package to use this function.", call. = FALSE)
  }
  checkmate::assert_numeric(log_p, any.missing = FALSE, len = length(log_q))
  checkmate::assert_numeric(log_q, any.missing = FALSE, len = length(log_p))
  checkmate::assert_matrix(log_liks, null.ok = TRUE, nrows = length(log_p))
  checkmate::assert_integerish(cores)
  checkmate::assert_flag(save_psis)

  approx_correction <- log_p - log_q

  if (is.null(log_liks)){
    # Handle underflow/overflow
    approx_correction <- approx_correction - max(approx_correction)
    log_ratios <- matrix(approx_correction, ncol = 1)
  } else {
    log_ratios <- -log_liks
    log_ratios <- log_ratios + approx_correction
    # Handle underflow/overflow
    log_ratio_max <- matrixStats::colMaxs(log_ratios)
    log_ratios <- sweep(log_ratios, MARGIN = 2, STATS = log_ratio_max)
  }
  psis_out <- psis.matrix(log_ratios, cores = cores, r_eff = rep(1, ncol(log_ratios)))

  if (is.null(log_liks)) {
    return(psis_out)
  }

  pointwise <- pointwise_loo_calcs(log_liks, psis_out)
  psis_loo_object(
    pointwise = pointwise,
    diagnostics = psis_out$diagnostics,
    dims = dim(psis_out),
    psis_object = if (save_psis) psis_out else NULL
  )
}

