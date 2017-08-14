# Effective sample size for PSIS
#
# @param w A vector or matrix (one column per observation) of normalized Pareto
#   smoothed weights (not log weights).
# @param rel_neff Relative ESS of exp(log_lik). rel_neff should be a scalar if w
#   is a vector and a vector of length ncol(w) if w is a matrix.
# @return A scalar if w is a vector. A vector of length ncol(w) if w is matrix.
#
psis_neff <- function(w, ...) {
  UseMethod("psis_neff")
}

psis_neff.default <- function(w, rel_neff = NULL, ...) {
  ss <- sum(w^2)
  if (is.null(rel_neff)) {
    warning("PSIS n_eff not adjusted based on MCMC n_eff.", call. = FALSE)
    return(1 / ss)
  }
  stopifnot(length(rel_neff) == 1)
  1 / ss * rel_neff
}

psis_neff.matrix <- function(w, rel_neff = NULL, ...) {
  ss <- colSums(w^2)
  if (is.null(rel_neff)) {
    warning("PSIS n_eff not adjusted based on MCMC n_eff.", call. = FALSE)
    return(1 / ss)
  }
  if (length(rel_neff) != length(ss))
    stop("rel_neff must have length ncol(w).", call. = FALSE)
  1 / ss * rel_neff
}


# Compute relative effective sample sizes
#
# @param x A vector, matrix, or 3-D array. In the case of the loo package, this
#   will be _exponentiated_ log-likelihood values.
# @param chain_id If x is not a 3-D array, a vector give chain indexes for each
#   value in x (if x is a vector) or each row of x (if x is a matrix).
# @return A scalar if x is a vector, or a vector if x is a matrix or 3-D array.
#
relative_neff <- function(x, ...) {
  UseMethod("relative_neff")
}

relative_neff.default <- function(x, chain_id, ...) {
  dim(x) <- c(length(x), 1)
  class(x) <- "matrix"
  relative_neff.matrix(x, chain_id)
}

relative_neff.matrix <- function(x, chain_id, ...) {
  x <- .llmat_to_array(x, chain_id)
  relative_neff.array(x)
}

relative_neff.array <- function(x, ...) {
  stopifnot(length(dim(x)) == 3)
  neff_vec <- apply(x, 3, .mcmc_neff)
  S <- prod(dim(x)[1:2]) # total samples = iter * chains
  neff_vec / S
}


# Convert (iter * chain) by obs matrix to iter by chain by obs array
#
# @param x matrix to convert.
# @param chain_id vector of chain ids.
# @return iter by chain by obs array
#
.llmat_to_array <- function(x, chain_id) {
  stopifnot(is.matrix(x))
  stopifnot(all(chain_id == as.integer(chain_id)))
  n_chain <- length(unique(chain_id))
  lldim <- dim(x)
  if (lldim[1] %% n_chain != 0) {
    stop("Number of rows in matrix not divisible ",
         "by number of chains in 'chain_id'.",
         call. = FALSE)
  }
  n_iter <- lldim[1] / n_chain
  n_obs <- lldim[2]

  a <- array(
    data = NA,
    dim = c(n_iter, n_chain, n_obs),
    dimnames = list(
      Iteration = NULL,
      Chain = NULL,
      Observation = NULL
    )
  )
  for (c in seq_len(n_chain)) {
    a[, c, ] <- x[chain_id == c, , drop = FALSE]
  }
  return(a)
}


# MCMC effective sample size calculation
#
# @param x An iterations by chains matrix of draws for a single parameter. In
#   the case of the loo package, this will be the _exponentiated_ log-likelihood
#   values for the ith observation.
# @return effective sample size based on rstan's calculation
#
.mcmc_neff <- function(x) {
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

  neff <- n_chain * n_iter
  if (rho_hat_sum > 0) {
    neff <- neff/(1 + 2 * rho_hat_sum)
  }

  return(neff)
}

# wrapper around stats::acf that returns only the info we need in .mcmc_neff
# @param x,lag_max Vector and integer passed to stats::acf
.acov <- function(x, lag_max) {
  cov <-
    stats::acf(x,
               lag.max = lag_max,
               plot = FALSE,
               type = "covariance")

  return(cov$acf[, , 1])
}


# Variance of the expected loo predictive density (EPD)
#
# Note: this is NOT for the expected loo _log_ predictive density (ELPD)
#
# @param w Vector or matrix of normalized Pareto smoothed weights
# @param w_un Vector of matrix of unnormalized Pareto smoothed weights
# @param rel_neff Precomputed relative effective sample size(s) of exp(log_lik).
# @return A scalar if w is a vector or a vector if w is a matrix.
#
var_epd <- function(w, w_un, rel_neff, ...) {
  UseMethod("var_epd")
}
var_epd.default <- function(w, w_un, rel_neff, ...) {
  stopifnot(length(w) == length(w_un),
            length(rel_neff) == 1)
  .var_epd_i(w, w_un, rel_neff)
}
var_epd.matrix <- function(w, w_un, rel_neff) {
  stopifnot(identical(dim(w), dim(w_un)),
            length(rel_neff) == ncol(w))
  sapply(1:ncol(w), function(i) {
    .var_epd_i(w[, i], w_un[, i], rel_neff[i])
  })
}

# @param w_i Vector of normalized weights for ith obs
# @param w_un_i Vector of unnormalized weights for ith obs
# @param rel_neff_i Scalar relative neff
# @param return A scalar
.var_epd_i <- function(w_i, w_un_i, rel_neff_i) {
  sum(w_i^2 * (1 / w_un_i - 1 / mean(w_un_i))^2) / rel_neff_i
}
