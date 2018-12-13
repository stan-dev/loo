#' Diagnostics for Laplace and ADVI approximations and Laplace-loo and ADVI-loo
#'
#' @param log_p The log-posterior (target) evaluated at samples from the proposal distribution (q).
#' @param log_q The log-density (proposal) evaluated at samples from the proposal distribution (q).
#' @param log_liks A log-likelihood matrix. See \code{\link{loo.matrix}}. Default is \code{NULL}. Then only the posterior is evaluated using the k_hat diagnostic.
#' @inheritParams loo
#' 
#' @return 
#' If log likelihoods are supplied, the function returns a \code{loo} object,
#' otherwise the function returns a \code{psis} object.
#' 
#' @seealso \code{loo} and \code{psis}
#'
#' @template loo-and-psis-references
#'
#' @export
psis_approximate_posterior <- function(log_p, log_q, log_liks = NULL, cores, save_psis, ...){
  checkmate::assert_numeric(log_p, any.missing = FALSE, len = length(log_q))
  checkmate::assert_numeric(log_q, any.missing = FALSE, len = length(log_p))
  checkmate::assert_matrix(log_liks, null.ok = TRUE, nrows = length(log_p))
  checkmate::assert_integerish(cores)
  checkmate::assert_flag(save_psis)

  approx_correction <- log_p - log_q
  approx_correction <- approx_correction - min(approx_correction) # Handle underflow/overflow
  
  if(is.null(log_liks)){
    log_ratios <- matrix(approx_correction, ncol = 1)
  } else {
    log_ratios <- -log_liks
    for(j in 1:ncol(log_ratios)){
      log_ratios[,j] <- log_ratios[,j] + approx_correction
    }
  }
  psis_out <- loo:::psis.matrix(log_ratios = log_ratios, cores = cores, r_eff = rep(1, ncol(log_ratios)))

  if(!is.null(log_liks)){
    pointwise <- loo:::pointwise_loo_calcs(log_liks, psis_out)
    obj <- loo:::psis_loo_object(pointwise = pointwise, diagnostics = psis_out$diagnostics, 
                                 dims = dim(psis_out), psis_object = if (save_psis) 
                                   psis_out else NULL)
    return(obj)
  } 
  return(psis_out)
}
