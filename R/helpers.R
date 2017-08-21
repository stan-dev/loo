#' @importFrom matrixStats logSumExp colLogSumExps colSums2 colVars

# more stable version of log(colMeans(exp(x)))
# @param x matrix
logColMeansExp <- function(x) {
  logS <- log(nrow(x))
  colLogSumExps(x) - logS
}

# Compute point estimates and standard errors from pointwise vectors
#
# @param x A matrix.
# @return An ncol(pointwise) by 2 matrix with columns 'Estimate' and 'SE'
#   and rownames equal to colnames(pointwise).
#
table_of_estimates <- function(x) {
  out <- cbind(
    Estimate = colSums2(x),
    SE = sqrt(nrow(x) * colVars(x))
  )
  rownames(out) <- colnames(x)
  return(out)
}


# Variance of the expected loo predictive density (EPD)
#
# Note: this is NOT for the expected loo _log_ predictive density (ELPD)
#
# @param w Vector or matrix of normalized Pareto smoothed weights
# @param w_un Vector of matrix of unnormalized Pareto smoothed weights
# @param rel_n_eff Precomputed relative effective sample size(s) of exp(log_lik).
# @return A scalar if w is a vector or a vector if w is a matrix.
#
var_epd <- function(w, w_un, rel_n_eff, ...) {
  UseMethod("var_epd")
}
var_epd.default <- function(w, w_un, rel_n_eff, ...) {
  stopifnot(length(w) == length(w_un),
            length(rel_n_eff) == 1)
  .var_epd_i(w, w_un, rel_n_eff)
}
var_epd.matrix <- function(w, w_un, rel_n_eff) {
  stopifnot(identical(dim(w), dim(w_un)),
            length(rel_n_eff) == ncol(w))
  sapply(1:ncol(w), function(i) {
    .var_epd_i(w[, i], w_un[, i], rel_n_eff[i])
  })
}

# @param w_i Vector of normalized weights for ith obs
# @param w_un_i Vector of unnormalized weights for ith obs
# @param rel_n_eff_i Scalar relative n_eff
# @param return A scalar
.var_epd_i <- function(w_i, w_un_i, rel_n_eff_i) {
  sum(w_i^2 * (1 / w_un_i - 1 / mean(w_un_i))^2) / rel_n_eff_i
}

# Effective sample size for PSIS
#
# @param w A vector or matrix (one column per observation) of normalized Pareto
#   smoothed weights (not log weights).
# @param rel_n_eff Relative ESS of exp(log_lik). rel_n_eff should be a scalar if
#   w is a vector and a vector of length ncol(w) if w is a matrix.
# @return A scalar if w is a vector. A vector of length ncol(w) if w is matrix.
#
psis_n_eff <- function(w, ...) {
  UseMethod("psis_n_eff")
}
psis_n_eff.default <- function(w, rel_n_eff = NULL, ...) {
  ss <- sum(w^2)
  if (is.null(rel_n_eff)) {
    warning("PSIS n_eff not adjusted based on MCMC n_eff.", call. = FALSE)
    return(1 / ss)
  }
  stopifnot(length(rel_n_eff) == 1)
  1 / ss * rel_n_eff
}
psis_n_eff.matrix <- function(w, rel_n_eff = NULL, ...) {
  ss <- colSums(w^2)
  if (is.null(rel_n_eff)) {
    warning("PSIS n_eff not adjusted based on MCMC n_eff.", call. = FALSE)
    return(1 / ss)
  }
  if (length(rel_n_eff) != length(ss))
    stop("rel_n_eff must have length ncol(w).", call. = FALSE)
  1 / ss * rel_n_eff
}


# Compute relative effective sample sizes
#
# @param x A vector, matrix, or 3-D array. In the case of the loo package, this
#   will be _exponentiated_ log-likelihood values.
# @param chain_id If x is not a 3-D array, a vector giving chain indexes for each
#   value in x (if x is a vector) or each row of x (if x is a matrix).
# @return A scalar if x is a vector, or a vector if x is a matrix or 3-D array.
#
relative_n_eff <- function(x, ...) {
  UseMethod("relative_n_eff")
}
relative_n_eff.default <- function(x, chain_id, ...) {
  dim(x) <- c(length(x), 1)
  class(x) <- "matrix"
  relative_n_eff.matrix(x, chain_id)
}
relative_n_eff.matrix <- function(x, chain_id, ...) {
  x <- llmatrix_to_array(x, chain_id)
  relative_n_eff.array(x)
}
relative_n_eff.array <- function(x, ...) {
  stopifnot(length(dim(x)) == 3)
  n_eff_vec <- apply(x, 3, mcmc_n_eff)
  S <- prod(dim(x)[1:2]) # total samples = iter * chains
  n_eff_vec / S
}


# MCMC effective sample size calculation
#
# @param x An iterations by chains matrix of draws for a single parameter. In
#   the case of the loo package, this will be the _exponentiated_ log-likelihood
#   values for the ith observation.
# @return MCMC effective sample size based on rstan's calculation.
#
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



# checking classes --------------------------------------------------------
is.psis <- function(x) {
  inherits(x, "psis") && is.list(x)
}
is.loo <- function(x) {
  inherits(x, "loo")
}
is.psis_loo <- function(x) {
  inherits(x, "psis_loo") && is.loo(x)
}
is.waic <- function(x) {
  inherits(x, "waic") && is.loo(x)
}


# validating and reshaping arrays/matrices  -------------------------------

# check for NAs and non-finite values in log-lik array/matrix/vector
# @param x log-lik array/matrix/vector
# @return x if no error is thrown.
#
validate_ll <- function(x) {
  stopifnot(!is.list(x), !anyNA(x), all(is.finite(x)))
  invisible(x)
}

# Convert iter by chain by obs array to (iter * chain) by obs matrix
#
# @param x array to convert.
# @return (iter * chain) by obs matrix
#
llarray_to_matrix <- function(x) {
  stopifnot(is.array(x), length(dim(x)) == 3)
  xdim <- dim(x)
  dim(x) <- c(prod(xdim[1:2]), xdim[3])
  unname(x)
}

# Convert (iter * chain) by obs matrix to iter by chain by obs array
#
# @param x matrix to convert.
# @param chain_id vector of chain ids.
# @return iter by chain by obs array
#
llmatrix_to_array <- function(x, chain_id) {
  stopifnot(is.matrix(x), all(chain_id == as.integer(chain_id)),
            length(chain_id) == nrow(x))

  lldim <- dim(x)
  chain_id <- as.integer(chain_id)
  n_chain <- length(unique(chain_id))
  if (max(chain_id) != n_chain) {
    stop("max(chain_id) not equal to the number of chains.",
         call. = FALSE)
  } else if (lldim[1] %% n_chain != 0) {
    stop("Number of rows in matrix not divisible ",
         "by number of chains in chain_id.",
         call. = FALSE)
  }

  n_iter <- lldim[1] / n_chain
  n_obs <- lldim[2]
  a <- array(data = NA, dim = c(n_iter, n_chain, n_obs))
  for (c in seq_len(n_chain)) {
    a[, c, ] <- x[chain_id == c, , drop = FALSE]
  }
  return(a)
}



#' Named lists
#'
#' Create a named list using specified names or, if names are omitted, using the
#' names of the objects in the list. The code \code{list(a = a, b = b)} becomes
#' \code{nlist(a,b)} and \code{list(a = a, b = 2)} becomes \code{nlist(a, b =
#' 2)}, etc.
#'
#' @export
#' @keywords internal
#' @param ... Objects to include in the list.
#' @return A named list.
#'
#' @seealso \code{\link[base]{list}}
#' @author Jonah Gabry
#' @examples
#'
#' # All variables already defined
#' a <- rnorm(100)
#' b <- mat.or.vec(10, 3)
#' nlist(a,b)
#'
#' # Define some variables in the call and take the rest from the environment
#' nlist(a, b, veggies = c("lettuce", "spinach"), fruits = c("banana", "papaya"))
#'
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }

  return(out)
}

# nocov start
# release reminders (for devtools)
release_questions <- function() {
  c(
    "Have you updated all references to the LOO paper?",
    "Have you updated inst/CITATION?",
    "Have you updated R code in vignette to match the code in the paper?"
  )
}
# nocov end
