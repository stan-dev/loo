#' Monte Carlo diagnostics of the CV bridge
#'
#' @export
#' @inheritParams brute_force_diags
#' @param log_lik_full An \eqn{S \times N} matrix of the log-likelihood at
#'   the full-data draws.
#' @param chain_id_full The chain of each full-data draw.
#' @return A `cv_diags` object. The total MCSE is the IJ over the shared
#'   full-data draws plus the variance of the refit-draw terms.
#'
#' @examples
#' set.seed(1)
#' y <- c(rnorm(19), 4)
#' N <- length(y)
#' S <- 1000
#' mu_full <- rnorm(S, mean(y), sqrt(1 / N))
#' log_lik_full <- sapply(y, function(y_i) dnorm(y_i, mu_full, log = TRUE))
#' log_lik_fold <- sapply(seq_len(N), function(i) {
#'   mu_i <- rnorm(S, mean(y[-i]), sqrt(1 / (N - 1)))
#'   dnorm(y[i], mu_i, log = TRUE)
#' })
#' cv_bridge_diags(log_lik_fold, log_lik_full)
#'
cv_bridge_diags <- function(x, log_lik_full, chain_id = NULL,
                            chain_id_full = chain_id) {
  .cv_bridge_diags(x, log_lik_full, chain_id = chain_id,
                   chain_id_full = chain_id_full, gap_na = 500)
}

#' CV bridge core with the gap cut-off
#'
#' @noRd
#' @param gap_na A fold with a cross-density gap above `gap_na` nat gets
#'   `NA`.
#'
.cv_bridge_diags <- function(x, log_lik_full, chain_id = NULL,
                             chain_id_full = chain_id, gap_na = 500) {
  ll <- .assert_log_lik(x)
  checkmate::assert_matrix(log_lik_full, mode = "numeric", nrows = nrow(ll),
                           ncols = ncol(ll))
  S <- nrow(ll)
  N <- ncol(ll)
  .assert_chain_id(chain_id, S)
  .assert_chain_id(chain_id_full, S)
  checkmate::assert_number(gap_na, lower = 0)

  folds <- lapply(seq_len(N), function(i) {
    .cv_bridge_i(ll[, i], log_lik_full[, i], chain_id, chain_id_full,
                 gap_na)
  })
  influence <- lapply(folds, `[[`, "influence")
  mcse <- vapply(folds, `[[`, numeric(1), "mcse")

  mcse_total <- NA_real_
  if (all(is.finite(mcse))) {
    psi_full <- Reduce(`+`, lapply(influence, `[[`, "full"))
    var_fold <- vapply(influence, function(x) .ij_variance(x$fold, chain_id), numeric(1))
    mcse_total <- sqrt(.ij_variance(psi_full, chain_id_full) + sum(var_fold))
  }

  pointwise <- data.frame(elpd = vapply(folds, `[[`, numeric(1), "elpd"),
                          mcse_elpd = mcse, bias_elpd = 0,
                          gap = vapply(folds, `[[`, numeric(1), "gap"),
                          converged = vapply(folds, `[[`, logical(1), "converged"),
                          stringsAsFactors = FALSE)
  .cv_diags(pointwise, mcse_total, list(influence = influence), "cv_bridge", S)
}

#' CV bridge estimate, MCSE and influences of one fold
#'
#' Draws theta_f^j (j = 1..S1) come from the full-data posterior, draws
#' theta_l^k (k = 1..S2) from the fit without the fold. The ratio of the two
#' unnormalized densities is rho(theta) = p(y_i | theta), so the inputs are
#' log-lik values only:
#'   log_rho_fold = log p(y_i | theta_l^k),  log_rho_full = log p(y_i | theta_f^j).
#' The fixed point r = N(r) / D(r) of the bridge estimator estimates
#' p(y_i | y_-i). The code holds D on the scale r * D.
#'
#' @noRd
#' @param log_rho_fold Log-lik of the fold at the refit draws.
#' @param log_rho_full Log-lik of the fold at the full-data draws.
#' @param chain_id,chain_id_full `NULL` or a valid `chain_id` of each set.
#' @param gap_na A gap above `gap_na` nat gives `NA`.
#' @param tol,max_iter Convergence of the fixed point on the log scale.
#' @return A list with `elpd`, `mcse`, `gap`, `converged` and `influence`
#'   (a list with `full` and `fold`, or `NULL`).
#'
.cv_bridge_i <- function(log_rho_fold, log_rho_full, chain_id, chain_id_full,
                         gap_na, tol = 1e-10, max_iter = 1000) {
  gap <- abs(stats::median(log_rho_fold) - stats::median(log_rho_full))
  failed <- list(elpd = NA_real_, mcse = NA_real_, gap = gap, converged = FALSE,
                 influence = NULL)
  if (!is.finite(gap) || gap > gap_na) {
    return(failed)
  }

  n_eff_full <- .n_eff_of_mean(log_rho_full, chain_id_full)
  n_eff_fold <- .n_eff_of_mean(log_rho_fold, chain_id)
  log_s1 <- log(n_eff_full / (n_eff_full + n_eff_fold))
  log_s2 <- log(n_eff_fold / (n_eff_full + n_eff_fold))

  log_r <- logMeanExp(log_rho_fold)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    terms <- .bridge_terms(log_rho_fold, log_rho_full, log_r, log_s1, log_s2)
    log_r_new <- log_r + logMeanExp(terms$log_n) - logMeanExp(terms$log_d)
    if (!is.finite(log_r_new)) {
      return(failed)
    }
    step <- abs(log_r_new - log_r)
    log_r <- log_r_new
    if (step < tol) {
      converged <- TRUE
      break
    }
  }

  terms <- .bridge_terms(log_rho_fold, log_rho_full, log_r, log_s1, log_s2)
  n <- exp(terms$log_n)                 # N_k, refit draws
  d <- exp(terms$log_d)                 # r * D_j, full-data draws
  d_fold <- exp(terms$log_d_fold)       # r * D at the refit draws
  s2 <- exp(log_s2)

  # per-fold MCSE: delta method on the ratio of two independent means
  rel_var <- .ij_variance((n - mean(n)) / (mean(n) * length(n)), chain_id) +
    .ij_variance((d - mean(d)) / (mean(d) * length(d)), chain_id_full)

  # influence of each draw on log r, from the implicit function theorem;
  # d_check includes the two feedback terms of the fixed point
  d_check <- mean(d) + s2 * mean(n * d_fold) - s2 * mean(d^2)
  influence <- list(
    full = -(d - mean(d)) / (d_check * length(d)),
    fold = (n - mean(n)) / (d_check * length(n))
  )

  list(elpd = log_r, mcse = sqrt(log1p(rel_var)), gap = gap, converged = converged,
       influence = influence)
}

#' Log summands of the bridge estimator at `log_r`
#'
#' @noRd
#' @return A list: `log_n` (refit draws), `log_d` (full-data draws, scaled
#'   by r) and `log_d_fold` (the scaled D summand at the refit draws).
#'
.bridge_terms <- function(log_rho_fold, log_rho_full, log_r, log_s1, log_s2) {
  log_den_fold <- matrixStats::colLogSumExps(rbind(log_s1 + log_rho_fold,
                                                    log_s2 + log_r))
  log_den_full <- matrixStats::colLogSumExps(rbind(log_s1 + log_rho_full,
                                                    log_s2 + log_r))
  list(
    log_n = log_rho_fold - log_den_fold,
    log_d = log_r - log_den_full,
    log_d_fold = log_r - log_den_fold
  )
}
