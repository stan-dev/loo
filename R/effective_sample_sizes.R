#' Convenience function for computing relative efficiencies
#'
#' \code{relative_eff} computes the the MCMC effective sample size divided by
#' the total sample size. The matrix and array methods perform the computations
#' independently for each observation.
#'
#' @export
#' @param x A vector, matrix, or 3-D array. See the \strong{Methods (by class)}
#'   section below for details on the shape of \code{x}. For use with the
#'   \code{loo} function, the values in \code{x} should be likelihood values
#'   (i.e., \code{exp(log_lik)}). For generic \code{use} with
#'   \code{\link{psis}}, the values in \code{x} should be the reciprocal of the
#'   importance ratios (i.e., \code{exp(-log_ratios)}).
#' @param chain_id A vector of length \code{NROW(x)} containing MCMC chain
#'   indexes for each each row of \code{x} (if a matrix) or each value in
#'   \code{x} (if a vector). No \code{chain_id} is needed if \code{x} is a 3-D
#'   array. If there are \code{C} chains then valid chain indexes are values
#'   in \code{1:C}.
#' @param ... Arguments passed to the methods.
#' @return A scalar if \code{x} is a vector, or a vector if \code{x} is a matrix
#'   or 3-D array.
#'
relative_eff <- function(x, ...) {
  UseMethod("relative_eff")
}

#' @export
#' @templateVar fn relative_eff
#' @template vector
#'
relative_eff.default <- function(x, chain_id, ...) {
  dim(x) <- c(length(x), 1)
  class(x) <- "matrix"
  relative_eff.matrix(x, chain_id)
}

#' @export
#' @templateVar fn relative_eff
#' @template matrix
#'
relative_eff.matrix <- function(x, chain_id, ...) {
  x <- llmatrix_to_array(x, chain_id)
  relative_eff.array(x)
}

#' @export
#' @templateVar fn relative_eff
#' @template array
#'
relative_eff.array <- function(x, ...) {
  stopifnot(length(dim(x)) == 3)
  n_eff_vec <- apply(x, 3, mcmc_n_eff)
  S <- prod(dim(x)[1:2]) # total samples = iter * chains
  n_eff_vec / S
}


#' Effective sample size for PSIS
#'
#' @noRd
#' @param w A vector or matrix (one column per observation) of normalized Pareto
#'   smoothed weights (not log weights).
#' @param r_eff Relative effective sample size of \code{exp(log_lik)} or
#'   \code{exp(-log_ratios)}. \cope{r_eff} should be a scalar if \code{w} is a
#'   vector and a vector of length \code{ncol(w)} if \code{w} is a matrix.
#' @return A scalar if \code{w} is a vector. A vector of length \code{ncol(w)}
#'   if \code{w} is matrix.
#'
psis_n_eff <- function(w, ...) {
  UseMethod("psis_n_eff")
}
psis_n_eff.default <- function(w, r_eff = NULL, ...) {
  ss <- sum(w^2)
  if (is.null(r_eff)) {
    warning("PSIS n_eff not adjusted based on MCMC n_eff.", call. = FALSE)
    return(1 / ss)
  }
  stopifnot(length(r_eff) == 1)
  1 / ss * r_eff
}
psis_n_eff.matrix <- function(w, r_eff = NULL, ...) {
  ss <- colSums(w^2)
  if (is.null(r_eff)) {
    warning("PSIS n_eff not adjusted based on MCMC n_eff.", call. = FALSE)
    return(1 / ss)
  }
  if (length(r_eff) != length(ss))
    stop("r_eff must have length ncol(w).", call. = FALSE)
  1 / ss * r_eff
}

#' MCMC effective sample size calculation
#'
#' @noRd
#' @param x An iterations by chains matrix of draws for a single parameter. In
#'   the case of the loo package, this will be the _exponentiated_ log-likelihood
#'   values for the ith observation.
#' @return MCMC effective sample size based on rstan's calculation.
#'
mcmc_n_eff <- function(x) {
  stopifnot(is.matrix(x))
  n_chain <- ncol(x)
  n_iter <- nrow(x)

  acov <- apply(x, 2, .acov, lag_max = n_iter - 1)
  chain_means <- colMeans(x)
  mean_var <- mean(acov[1, ]) * n_iter / (n_iter - 1)
  var_plus <- mean_var * (n_iter - 1) / n_iter
  if (n_chain > 1) {
    var_plus <- var_plus + var(chain_means)
  }

  rho_hat_sum <- 0
  for (t in 2:nrow(acov)) {
    rho_hat <- 1 - (mean_var - mean(acov[t, ])) / var_plus
    if (is.nan(rho_hat))
      rho_hat <- 0
    if (rho_hat < 0)
      break
    rho_hat_sum <- rho_hat_sum + rho_hat
  }

  n_eff <- n_chain * n_iter
  if (rho_hat_sum > 0)
    n_eff <- n_eff / (1 + 2 * rho_hat_sum)

  return(n_eff)
}


# wrapper around stats::acf that returns only the info we need in mcmc_n_eff
# @param x,lag_max Vector and integer passed to stats::acf
.acov <- function(x, lag_max) {
  cov <-
    stats::acf(x,
               lag.max = lag_max,
               plot = FALSE,
               type = "covariance")

  return(cov$acf[, , 1])
}


