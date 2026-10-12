#' Monte Carlo diagnostics of brute-force LOO
#'
#' @export
#' @param x An \eqn{S \times N} matrix. Column \eqn{i} is the
#'   log-likelihood of fold \eqn{i} at the draws of the fit without fold
#'   \eqn{i}.
#' @param log_lik_full An optional \eqn{S \times N} matrix of the
#'   log-likelihood at the full-data draws. Its largest density per fold
#'   bounds the tail of the truncated-GPD MCSE.
#' @param chain_id The chain of each refit draw. `NULL`: the draws are independent.
#' @return A `cv_diags` object. The total MCSE is the root sum of squares:
#'   the refits have independent draws.
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
#' brute_force_diags(log_lik_fold, log_lik_full)
#'
brute_force_diags <- function(x, log_lik_full = NULL, chain_id = NULL) {
  .brute_force_diags(x, log_lik_full = log_lik_full,
                     chain_id = chain_id, upper = Inf)
}

#' Brute-force core with the family bound
#'
#' @noRd
#' @param upper Bound of \eqn{p(y_i \mid \theta)} per fold, or `Inf`.
#'
.brute_force_diags <- function(x, log_lik_full = NULL, chain_id = NULL,
                               upper = Inf) {
  ll <- .assert_log_lik(x)
  S <- nrow(ll)
  N <- ncol(ll)
  .assert_chain_id(chain_id, S)
  if (!is.null(log_lik_full)) {
    checkmate::assert_matrix(log_lik_full, mode = "numeric", nrows = S,
                             ncols = N)
  }
  checkmate::assert_numeric(upper, lower = 0, any.missing = FALSE)
  upper <- rep_len(upper, N)
  if (!is.null(log_lik_full)) {
    upper <- pmin(upper, exp(matrixStats::colMaxs(log_lik_full)))
  }

  ll_array <- if (!is.null(chain_id)) llmatrix_to_array(ll, chain_id)
  per_fold <- lapply(seq_len(N), function(i) {
    ll_i <- ll[, i]
    if (!is.null(chain_id)) {
      ll_i <- matrix(ll_array[, , i], nrow = dim(ll_array)[1])
    }
    .log_mean_diags(ll_i, are_log = TRUE, upper = upper[i])
  })
  mcse <- vapply(per_fold, `[[`, numeric(1), "mcse")
  pointwise <- data.frame(elpd = colLogMeanExps(ll),
                          mcse_elpd = mcse,
                          bias_elpd = vapply(per_fold, `[[`, numeric(1), "bias"),
                          stringsAsFactors = FALSE)
  .cv_diags(pointwise, sqrt(sum(mcse^2)), list(), "brute_force", S)
}
