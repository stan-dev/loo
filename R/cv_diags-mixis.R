#' Monte Carlo diagnostics of mixture importance sampling (MixIS)
#'
#' The draws come from the mixture \eqn{q(\theta) \propto p(\theta \mid y)
#' \sum_j \alpha_j / p(y_j \mid \theta)}. With \eqn{m_s = \sum_j \alpha_j /
#' p(y_j \mid \theta^s)}, the weights are \eqn{u_s = 1 / m_s} for the
#' posterior and \eqn{v_{i,s} = u_s / p(y_i \mid \theta^s)} for fold \eqn{i}.
#' The estimate is \eqn{\widehat{\mathrm{elpd}}_i = \log \bar u - \log \bar v_i}.
#'
#' @export
#' @param x An \eqn{S \times N} matrix of the log-likelihood at the
#'   mixture draws.
#' @param alpha The mixture weights that sampled `x`. `NULL`:
#'   \eqn{\alpha_j = 1}.
#' @param chain_id The chain of each draw. `NULL`: the draws are independent.
#' @return A `cv_diags` object. The total MCSE is the joint total of
#'   \eqn{u} and all \eqn{v_i}.
#'
#' @examples
#' # real use needs draws from the mixture fitted with `alpha`
#' LLmat <- example_loglik_matrix()
#' chain_id <- rep(1:2, each = nrow(LLmat) / 2)
#' mixis_diags(LLmat, chain_id = chain_id)
#'
mixis_diags <- function(x, alpha = NULL, chain_id = NULL) {
  ll <- .assert_log_lik(x)
  S <- nrow(ll)
  N <- ncol(ll)
  .assert_chain_id(chain_id, S)
  alpha <- .validate_alpha(alpha, N)

  log_u <- -.mixis_log_m(ll, alpha)
  log_v <- log_u - ll
  elpd <- logMeanExp(log_u) - colLogMeanExps(log_v)

  u <- exp(log_u - max(log_u))
  v <- exp(sweep(log_v, 2, matrixStats::colMaxs(log_v)))
  rel_u <- (u - mean(u)) / mean(u)
  rel_v <- sweep(sweep(v, 2, colMeans(v)), 2, colMeans(v), "/")

  var_u <- .ij_variance(rel_u / S, chain_id)
  var_v <- apply(rel_v, 2, function(x) .ij_variance(x / S, chain_id))
  mcse <- apply(rel_u - rel_v, 2, function(x) sqrt(.ij_variance(x / S, chain_id)))
  mcse_total <- sqrt(.ij_variance((N * rel_u - rowSums(rel_v)) / S, chain_id))

  pointwise <- data.frame(
    elpd = elpd, mcse_elpd = mcse,
    bias_elpd = var_v / 2 - var_u / 2,
    pareto_k_v = apply(log_v, 2, .right_tail_khat),
    stringsAsFactors = FALSE
  )
  diagnostics <- list(
    pareto_k_u = .right_tail_khat(log_u),
    n_eff_u = psis_n_eff(exp(log_u - matrixStats::logSumExp(log_u)))
  )
  .cv_diags(pointwise, mcse_total, diagnostics, "mixis", S)
}

#' Balanced mixture weights from a PSIS result
#'
#' \eqn{\alpha_j = \exp(e_j - \max_k e_k)} with the PSIS elpd \eqn{e_j}.
#'
#' @export
#' @param x A [psis_diags()] result, or a [loo()] result.
#' @return A vector of \eqn{N} weights; the largest is 1.
#'
#' @examples
#' LLmat <- example_loglik_matrix()
#' summary(mixis_alpha(psis_diags(LLmat)))
#'
mixis_alpha <- function(x) {
  elpd <- if (inherits(x, "cv_diags") && identical(x$method, "psis")) {
    x$pointwise$elpd
  } else if (inherits(x, "psis_loo")) {
    x$pointwise[, "elpd_loo"]
  } else {
    stop("'x' must be a psis_diags() or a loo() result.", call. = FALSE)
  }
  exp(elpd - max(elpd))
}

#' Pre-flight check: the predicted ESS of the MixIS weights
#'
#' The relative ESS of importance sampling from the posterior to the
#' mixture is \eqn{1 / (E[m_\alpha]\, E[1 / m_\alpha])}, with both
#' expectations over the posterior draws.
#'
#' @export
#' @param x An \eqn{S \times N} matrix of the log-likelihood at the
#'   posterior draws.
#' @param alpha The candidate mixture weights.
#' @param ndraws The number of mixture draws to sample. `NULL`: \eqn{S}.
#' @return A list: `r_eff`, `n_eff_u`, `pareto_k` (above 0.7 the prediction is
#'   not reliable) and `ratio` (`r_eff` of `alpha` over that of
#'   \eqn{\alpha_j = 1}).
#'
#' @examples
#' LLmat <- example_loglik_matrix()
#' mixis_n_eff(LLmat, alpha = mixis_alpha(psis_diags(LLmat)))
#'
mixis_n_eff <- function(x, alpha, ndraws = NULL) {
  .assert_log_lik(x)
  checkmate::assert_count(ndraws, positive = TRUE, null.ok = TRUE)
  N <- ncol(x)
  alpha <- .validate_alpha(alpha, N)
  if (is.null(ndraws)) {
    ndraws <- nrow(x)
  }
  balanced <- .predicted_mixis_n_eff(x, alpha, ndraws)
  plain <- .predicted_mixis_n_eff(x, rep(1, N), ndraws)
  c(balanced, list(ratio = balanced$r_eff / plain$r_eff))
}

#' \eqn{\log m_s = \log \sum_j \alpha_j / p(y_j \mid \theta^s)}
#' @noRd
.mixis_log_m <- function(ll, alpha) {
  matrixStats::rowLogSumExps(sweep(-ll, 2, log(alpha), "+"))
}

#' @noRd
.predicted_mixis_n_eff <- function(ll, alpha, ndraws) {
  log_m <- .mixis_log_m(ll, alpha)
  r_eff <- exp(-(logMeanExp(log_m) + logMeanExp(-log_m)))
  list(r_eff = r_eff, n_eff_u = ndraws * r_eff, pareto_k = .right_tail_khat(log_m))
}

#' Check the mixture weights; `NULL` gives \eqn{\alpha_j = 1}
#' @noRd
.validate_alpha <- function(alpha, N) {
  if (is.null(alpha)) {
    return(rep(1, N))
  }
  checkmate::assert_numeric(alpha, lower = 0, finite = TRUE, any.missing = FALSE,
                            len = N)
  if (!any(alpha > 0)) {
    stop("At least one value of 'alpha' must be positive.", call. = FALSE)
  }
  alpha
}
