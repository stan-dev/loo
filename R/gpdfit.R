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
#' @return A named list with components `k_hat`, `sigma_hat`, `k`, `k_w`, and `k_d`
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
  # see section 4 of Zhang and Stephens (2009)
  if (sort_x) {
    x <- sort.int(x)
  }
  N <- length(x)
  prior <- 3
  M <- min_grid_pts + floor(sqrt(N))
  jj <- seq_len(M)
  xstar <- x[floor(N / 4 + 0.5)] # first quartile of sample
  theta <- 1 / x[N] + (1 - sqrt(M / (jj - 0.5))) / prior / xstar
  k <- matrixStats::rowMeans2(log1p(-theta %o% x))
  l_theta <- N * (log(-theta / k) - k - 1) # profile log-lik
  w_theta <- exp(l_theta - matrixStats::logSumExp(l_theta)) # normalize
  theta_hat <- sum(theta * w_theta)
  k_hat <- mean.default(log1p(-theta_hat * x))
  sigma_hat <- -k_hat / theta_hat

  # quadrature weights for k are same as for theta
  k_w <- w_theta
  # quadrature weights are just the normalized likelihoods
  # we get the unnormalized posterior by multiplying these by the prior
  k_d <- k_w * dgpd(-theta, mu = -1 / x[N], sigma = 1 / prior / xstar, k = 0.5)
  # normalize using the trapezoidal rule
  Z <- sum((k_d[-M] + k_d[-1]) * (k[-M] - k[-1])) / 2
  k_d <- k_d / Z

  # adjust k_hat based on weakly informative prior, Gaussian centered on 0.5.
  # this stabilizes estimates for very small Monte Carlo sample sizes and low neff
  if (wip) {
    k_hat <- (k_hat * N + 0.5 * 10) / (N + 10)
  }

  if (is.na(k_hat)) {
    k_hat <- Inf
    sigma_hat <- NaN
  }

  nlist(k_hat, sigma_hat, k, k_w, k_d)
}
