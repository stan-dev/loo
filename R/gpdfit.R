#' Estimate parameters of the Generalized Pareto distribution
#'
#' Estimate the parameters \eqn{k} and \eqn{\sigma} of the generalized Pareto
#' distribution (assuming location parameter is 0), given a sample \eqn{x}. By
#' default the Pareto fit uses a prior for \eqn{k}, which will stabilize
#' estimates for very small sample sizes and low effective sample sizes in the
#' case of MCMC samples. The weakly informative prior is a Gaussian prior
#' centered at 0.5.
#'
#' @export
#' @param x A numeric vector. The sample from which to estimate the parameters.
#' @param wip Logical indicating whether to adjust \eqn{k} based on a weakly
#'   informative Gaussian prior centered on 0.5. Defaults to \code{TRUE}.
#' @param min_grid_pts The minimum number of grid points used in the fitting
#'   algorithm. The actual number used is \code{min_grid_pts +
#'   floor(sqrt(length(x)))}.
#' @param sort_x If \code{TRUE} (the default), the first step in the fitting
#'   algorithm is to sort the elements of \code{x}. If \code{x} is already
#'   sorted in ascending order then \code{sort_x} can be set to \code{FALSE} to
#'   skip the initial sorting step.
#' @return A named list with components \code{k} and \code{sigma}.
#'
#' @details Here the parameter \eqn{k} is the negative of \eqn{k} in Zhang &
#'   Stephens (2009).
#'
#' @seealso \code{\link{psis}}, \code{\link{pareto-k-diagnostic}},
#'   \code{\link{loo-package}}
#'
#' @references
#' Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
#' for the generalized Pareto distribution. \emph{Technometrics} \strong{51},
#' 316-325.
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
  w_theta <- 1 / vapply(jj, FUN.VALUE = 0, FUN = function(j) {
    sum(exp(l_theta - l_theta[j]))
  })
  theta_hat <- sum(theta * w_theta)
  k <- mean.default(log1p(-theta_hat * x))
  sigma <- -k / theta_hat

  if (wip) {
    k <- adjust_k_wip(k, n = N)
  }

  if (is.nan(k)) {
    k <- Inf
  }

  nlist(k, sigma)
}


# internal ----------------------------------------------------------------

lx <- function(a,x) {
  a <- -a
  k <- sapply(a, FUN = function(y) mean(log1p(y * x)))
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
  nplusa <- n + a
  k * n / nplusa + a * 0.5 / nplusa
}


#' Inverse CDF of generalized pareto distribution
#' (assuming location parameter is 0)
#'
#' @noRd
#' @param p Vector of probabilities.
#' @param k Scalar shape parameter.
#' @param sigma Scalar scale parameter.
#' @return Vector of quantiles.
#'
qgpd <- function(p, k, sigma) {
  if (is.nan(sigma) || sigma <= 0)
    return(rep(NaN, length(p)))

  sigma * expm1(-k * log1p(-p)) / k
}
