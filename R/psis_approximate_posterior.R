#' Diagnostics for Laplace and ADVI approximations and Laplace-loo and ADVI-loo
#'
#' @param log_p The log-posterior (target) evaluated at S samples from the
#'   proposal distribution (g). A vector of length S.
#' @param log_g The log-density (proposal) evaluated at S samples from the
#'   proposal distribution (g). A vector of length S.
#' @param log_q Deprecated argument name (the same as log_g).
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
psis_approximate_posterior <- function(log_p = NULL, log_g = NULL, log_liks = NULL,
                                       cores, save_psis, ..., log_q = NULL) {
  if (!is.null(log_q)) {
    .Deprecated(msg = "psis_approximate_posterior() argument log_q has been changed to log_g")
    log_g <- log_q
  }
  checkmate::assert_numeric(log_p, any.missing = FALSE, len = length(log_g), null.ok = FALSE)
  checkmate::assert_numeric(log_g, any.missing = FALSE, len = length(log_p), null.ok = FALSE)
  checkmate::assert_matrix(log_liks, null.ok = TRUE, nrows = length(log_p))
  checkmate::assert_integerish(cores)
  checkmate::assert_flag(save_psis)

  if (is.null(log_liks)) {
    approx_correction <- log_p - log_g
    # Handle underflow/overflow
    approx_correction <- approx_correction - max(approx_correction)
    log_ratios <- matrix(approx_correction, ncol = 1)
  } else {
    log_ratios <- correct_log_ratios(log_ratios = -log_liks, log_p = log_p, log_g = log_g)
  }
  psis_out <- psis.matrix(log_ratios, cores = cores, r_eff = rep(1, ncol(log_ratios)))

  if (is.null(log_liks)) {
    return(psis_out)
  }

  pointwise <- pointwise_loo_calcs(log_liks, psis_out)
  importance_sampling_loo_object(
    pointwise = pointwise,
    diagnostics = psis_out$diagnostics,
    dims = dim(psis_out),
    is_method = "psis",
    is_object = if (save_psis) psis_out else NULL
  )
}



#' Correct log ratios for posterior approximations
#'
#' @inheritParams psis_approximate_posterior
#' @inheritParams ap_psis
#' @noRd
#' @keywords internal
correct_log_ratios <- function(log_ratios, log_p, log_g) {
  approx_correction <- log_p - log_g
  log_ratios <- log_ratios + approx_correction
  # Handle underflow/overflow
  log_ratio_max <- apply(log_ratios, 2, max)
  log_ratios <- sweep(log_ratios, MARGIN = 2, STATS = log_ratio_max)

  log_ratios
}


#' Pareto smoothed importance sampling (PSIS)
#' using approximate posteriors
#' @inheritParams psis_approximate_posterior
#' @param log_ratios The log-likelihood ratios (ie -log_liks)
#' @param ... Currently not in use.
ap_psis <- function(log_ratios, log_p, log_g, ...) {
  UseMethod("ap_psis")
}

#' @export
#' @templateVar fn ap_psis
#' @template array
#'
ap_psis.array <-
  function(log_ratios, log_p, log_g, ...,
           cores = getOption("mc.cores", 1)) {
    cores <- loo_cores(cores)
    stopifnot(length(dim(log_ratios)) == 3)
    log_ratios <- validate_ll(log_ratios)
    log_ratios <- llarray_to_matrix(log_ratios)

    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    ap_psis.matrix(log_ratios = log_ratios,
                   log_p = log_p,
                   log_g = log_g,
                   cores = 1)
  }

#' @export
#' @templateVar fn ap_psis
#' @template matrix
#'
ap_psis.matrix <- function(log_ratios, log_p, log_g,
           ...,
           cores = getOption("mc.cores", 1)) {
    checkmate::assert_numeric(log_p, len = nrow(log_ratios))
    checkmate::assert_numeric(log_g, len = nrow(log_ratios))
    cores <- loo_cores(cores)
    log_ratios <- validate_ll(log_ratios)

    log_ratios <- correct_log_ratios(log_ratios, log_p = log_p, log_g = log_g)

    do_psis(log_ratios, r_eff = rep(1, ncol(log_ratios)), cores = cores)
  }

#' @export
#' @templateVar fn ap_psis
#' @template vector
#'
ap_psis.default <- function(log_ratios, log_p, log_g, ...) {
    stopifnot(is.null(dim(log_ratios)) || length(dim(log_ratios)) == 1)
    dim(log_ratios) <- c(length(log_ratios), 1)
    warning("llfun values do not return a matrix, coerce to matrix")
    ap_psis.matrix(as.matrix(log_ratios), log_p, log_g, cores = 1)
  }


