#' Monte Carlo diagnostics of PSIS-LOO
#'
#' @export
#' @param x An \eqn{S \times N} matrix of the log-likelihood at the
#'   posterior draws. A column is a fold: one observation (LOO) or the sum
#'   over a group (LOGO).
#' @param chain_id The chain of each draw. `NULL`: the draws are independent.
#' @param r_eff Relative efficiency of each fold, as from [relative_eff()].
#'   `NULL`: computed from `chain_id`, or 1 without `chain_id`.
#' @param psis_object An existing [psis()] result for `-x`.
#' @return A `cv_diags` object. The per-fold MCSE is the truncated-GPD MCSE
#'   of the raw ratios; the total MCSE is the infinitesimal jackknife (IJ).
#'
#' @examples
#' LLmat <- example_loglik_matrix()
#' chain_id <- rep(1:2, each = nrow(LLmat) / 2)
#' psis_diags(LLmat, chain_id = chain_id)
#'
psis_diags <- function(x, chain_id = NULL, r_eff = NULL, psis_object = NULL) {
  ll <- .assert_log_lik(x)
  S <- nrow(ll)
  N <- ncol(ll)
  .assert_chain_id(chain_id, S)
  if (is.null(r_eff)) {
    r_eff <- if (is.null(chain_id)) {
      1
    } else {
      likelihood <- exp(sweep(ll, 2, matrixStats::colMaxs(ll)))
      relative_eff(likelihood, chain_id = chain_id)
    }
  }
  checkmate::assert_numeric(r_eff, lower = 0, finite = TRUE, any.missing = FALSE)
  if (!length(r_eff) %in% c(1, N)) {
    stop("'r_eff' must have length 1 or ncol(x).", call. = FALSE)
  }
  r_eff <- rep_len(r_eff, N)
  if (is.null(psis_object)) {
    psis_object <- psis(-ll, r_eff = r_eff)
  } else if (!identical(as.integer(dim(psis_object)), as.integer(dim(ll)))) {
    stop("'psis_object' must have the dimensions of 'x'.", call. = FALSE)
  }

  lw <- weights(psis_object, log = TRUE)
  elpd <- matrixStats::colLogSumExps(ll + lw)

  per_fold <- lapply(seq_len(N), function(i) {
    .log_mean_diags(-ll[, i], are_log = TRUE, r_eff = r_eff[i])
  })
  mcse <- vapply(per_fold, `[[`, numeric(1), "mcse")

  # influence of draw s on elpd_i: w_is * (p_is / p_i - 1), summed over folds
  influence <- exp(sweep(ll + lw, 2, elpd)) - exp(lw)
  mcse_total <- sqrt(.ij_variance(rowSums(influence), chain_id))

  pareto_k <- psis_object$diagnostics$pareto_k
  pointwise <- data.frame(elpd = elpd, mcse_elpd = mcse,
                          bias_elpd = mcse^2 / 2, pareto_k = pareto_k,
                          stringsAsFactors = FALSE)
  diagnostics <- list(pareto_k = pareto_k, n_eff = psis_object$diagnostics$n_eff,
                      r_eff = r_eff)
  .cv_diags(pointwise, mcse_total, diagnostics, "psis", S)
}
