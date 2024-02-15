#' Estimate parameters of the Generalized Pareto distribution
#'
#' Given a sample \eqn{x}, Estimate the parameters \eqn{k} and \eqn{\sigma} of
#' the generalized Pareto distribution (GPD), assuming the location parameter is
#' 0. By default the fit uses a prior for \eqn{k}, which will stabilize
#' estimates for very small sample sizes (and low effective sample sizes in the
#' case of MCMC samples). The weakly informative prior is a Gaussian prior
#' centered at 0.5.
#'
#' @export
#' @param x A numeric vector. The sample from which to estimate the parameters.
#' @param wip Logical indicating whether to adjust \eqn{k} based on a weakly
#'   informative Gaussian prior centered on 0.5. Defaults to `TRUE`.
#' @param min_grid_pts The minimum number of grid points used in the fitting
#'   algorithm. The actual number used is `min_grid_pts + floor(sqrt(length(x)))`.
#' @param sort_x If `TRUE` (the default), the first step in the fitting
#'   algorithm is to sort the elements of `x`. If `x` is already
#'   sorted in ascending order then `sort_x` can be set to `FALSE` to
#'   skip the initial sorting step.
#' @return A named list with components `k` and `sigma`.
#'
#' @details Here the parameter \eqn{k} is the negative of \eqn{k} in Zhang &
#'   Stephens (2009).
#'
#' @seealso [psis()], [pareto-k-diagnostic]
#'
#' @references
#' Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
#' for the generalized Pareto distribution. *Technometrics* **51**, 316-325.
#'
gpdfit <- function(x, wip = TRUE, min_grid_pts = 30, sort_x = TRUE) {
  # See section 4 of Zhang and Stephens (2009)
  if (sort_x) {
    x <- sort.int(x)
  }
  N <- length(x)
  prior <- 3
  M <- min_grid_pts + floor(sqrt(N))
  jj <- seq_len(M)
  xstar <- x[floor(N / 4 + 0.5)] # first quartile of sample
  theta <- 1 / x[N] + (1 - sqrt(M / (jj - 0.5))) / prior / xstar
  l_theta <- N * lx(theta, x) # profile log-lik
  w_theta <- exp(l_theta - matrixStats::logSumExp(l_theta)) # normalize
  theta_hat <- sum(theta * w_theta)
  k <- mean.default(log1p(-theta_hat * x))
  sigma <- -k / theta_hat

  if (wip) {
    k <- adjust_k_wip(k, n = N)
  }

  if (is.na(k)) {
    k <- Inf
  }

  nlist(k, sigma)
}


# internal ----------------------------------------------------------------

lx <- function(a,x) {
  a <- -a
  k <- vapply(a, FUN = function(a_i) mean(log1p(a_i * x)), FUN.VALUE = numeric(1))
  log(a / k) - k - 1
}

#' Adjust k based on weakly informative prior, Gaussian centered on 0.5. This
#' will stabilize estimates for very small Monte Carlo sample sizes and low neff
#' cases.
#'
#' @noRd
#' @param k Scalar khat estimate.
#' @param n Integer number of tail samples used to fit GPD.
#' @return Scalar adjusted khat estimate.
#'
adjust_k_wip <- function(k, n) {
  a <- 10
  n_plus_a <- n + a
  k * n / n_plus_a + a * 0.5 / n_plus_a
}


#' Inverse CDF of generalized Pareto distribution
#' (assuming location parameter is 0)
#'
#' @noRd
#' @param p Vector of probabilities.
#' @param k Scalar shape parameter.
#' @param sigma Scalar scale parameter.
#' @return Vector of quantiles.
#'
qgpd <- function(p, k, sigma) {
  if (is.nan(sigma) || sigma <= 0) {
    return(rep(NaN, length(p)))
  }

  sigma * expm1(-k * log1p(-p)) / k
}
