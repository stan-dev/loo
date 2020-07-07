#' Split moment matching for efficient approximate leave-one-out cross-validation (LOO)
#'
#' A function that computes the split moment matching importance sampling loo.
#' Takes in the moment matching total transformation, transforms only half
#' of the draws, and computes a single elpd using multiple importance sampling.
#'
#' @param x A fitted model object.
#' @param upars A matrix containing the model parameters in unconstrained space
#'   where they can have any real value.
#' @param cov Logical; Indicate whether to match the covariance matrix of the
#'   samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched.
#' @param total_shift A vector representing the total shift made by the moment
#'   matching algorithm.
#' @param total_scaling A vector representing the total scaling of marginal
#'   variance made by the moment matching algorithm.
#' @param total_mapping A vector representing the total covariance
#'   transformation made by the moment matching algorithm.
#' @param i Observation index.
#' @param log_prob_upars A function that takes arguments `x` and `upars` and
#'   returns a matrix of log-posterior density values of the unconstrained
#'   posterior draws passed via `upars`.
#' @param log_lik_i_upars A function that takes arguments `x`, `upars`, and `i`
#'   and returns a vector of log-likeliood draws of the `i`th observation based
#'   on the unconstrained posterior draws passed via `upars`.
#' @param r_eff_i MCMC relative effective sample size of the `i`'th log
#'   likelihood draws.
#' @template cores
#' @template is_method
#' @param ... Further arguments passed to the custom functions documented above.
#'
#' @return A list containing the updated log-importance weights and
#' log-likelihood values. Also returns the updated MCMC effective sample size
#' and the integrand-specific log-importance weights.
#'
#'
#' @seealso [loo()], [loo_moment_match()]
#' @template moment-matching-references
#'
#'
loo_moment_match_split <- function(x, upars, cov, total_shift, total_scaling,
                     total_mapping, i, log_prob_upars,
                     log_lik_i_upars, r_eff_i, cores,
                     is_method, ...) {
  S <- dim(upars)[1]
  S_half <- as.integer(0.5 * S)
  mean_original <- colMeans(upars)

  # accumulated affine transformation
  upars_trans <- sweep(upars, 2, mean_original, "-")
  upars_trans <- sweep(upars_trans, 2, total_scaling, "*")
  if (cov) {
    upars_trans <- tcrossprod(upars_trans, total_mapping)
  }
  upars_trans <- sweep(upars_trans, 2, total_shift + mean_original, "+")
  attributes(upars_trans) <- attributes(upars)

  # inverse accumulated affine transformation
  upars_trans_inv <- sweep(upars, 2, mean_original, "-")
  if (cov) {
    upars_trans_inv <- tcrossprod(upars_trans_inv, solve(total_mapping))
  }
  upars_trans_inv <- sweep(upars_trans_inv, 2, total_scaling, "/")
  upars_trans_inv <- sweep(upars_trans_inv, 2, mean_original - total_shift, "+")
  attributes(upars_trans_inv) <- attributes(upars)

  # first half of upars_trans_half are T(theta)
  # second half are theta
  upars_trans_half <- upars
  take <- seq_len(S_half)
  upars_trans_half[take, ] <- upars_trans[take, , drop = FALSE]

  # first half of upars_half_inv are theta
  # second half are T^-1 (theta)
  upars_trans_half_inv <- upars
  take <- seq_len(S)[-seq_len(S_half)]
  upars_trans_half_inv[take, ] <- upars_trans_inv[take, , drop = FALSE]

  # compute log likelihoods and log probabilities
  log_prob_half_trans <- log_prob_upars(x, upars = upars_trans_half, ...)
  log_prob_half_trans_inv <- log_prob_upars(x, upars = upars_trans_half_inv, ...)
  log_liki_half <- log_lik_i_upars(x, upars = upars_trans_half, i = i, ...)

  # compute weights
  lwi_half <- -log_liki_half + log_prob_half_trans -
    (log_prob_half_trans +
       log(1 + exp(log_prob_half_trans_inv - log(prod(total_scaling)) -
                     log(det(total_mapping)) - log_prob_half_trans)))

  is_obj_half <- suppressWarnings(importance_sampling.default(lwi_half,
                                                              method = is_method,
                                                              r_eff = r_eff_i,
                                                              cores = cores))
  lwi_half <- as.vector(weights(is_obj_half))

  is_obj_f_half <- suppressWarnings(importance_sampling.default(lwi_half +
                                                                  log_liki_half,
                                                                method = is_method,
                                                                r_eff = r_eff_i,
                                                                cores = cores))
  lwfi_half <- as.vector(weights(is_obj_f_half))

  # relative_eff recomputation
  # currently ignores chain information
  # since we have two proposal distributions
  # compute S_eff separately from both and take the smaller
  take <- seq_len(S)[-seq_len(S_half)]
  log_liki_half_1 <- log_liki_half[take, drop = FALSE]
  dim(log_liki_half_1) <- c(length(take), 1, 1)
  take <- seq_len(S)[-seq_len(S_half)]
  log_liki_half_2 <- log_liki_half[take, drop = FALSE]
  dim(log_liki_half_2) <- c(length(take), 1, 1)
  r_eff_i1 <- loo::relative_eff(exp(log_liki_half_1), cores = cores)
  r_eff_i2 <- loo::relative_eff(exp(log_liki_half_2), cores = cores)
  r_eff_i <- min(r_eff_i1,r_eff_i2)

  list(
    lwi = lwi_half,
    lwfi = lwfi_half,
    log_liki = log_liki_half,
    r_eff_i = r_eff_i
  )
}
