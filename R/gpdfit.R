#' Estimate parameters of the Generalized Pareto distribution
#'
#' Estimate the parameters \eqn{k} and \eqn{\sigma} of the generalized Pareto
#' distribution (assuming location parameter is 0), given a sample \eqn{x}. The
#' Pareto fit uses a prior for \eqn{k}, which will stabilize estimates for very
#' small Monte Carlo sample sizes and low effective sample sizes. The weakly
#' informative prior is a Gaussian prior centered on 0.5.
#'
#' @export
#' @param x A numeric vector. The sample from which to estimate the parameters.
#' @param wip Logical indicating whether to adjust \eqn{k} based on a weakly
#'   informative Gaussian prior centered on 0.5. Defaults to \code{TRUE}.
#' @param min_grid_pts The minimum number of grid points used. The actual number
#'   used is \code{min_grid_pts + floor(sqrt(length(x)))}.
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
gpdfit <- function(x, wip = TRUE, min_grid_pts = 30) {
  N <- length(x)
  x <- sort.int(x, method = "quick")
  prior <- 3
  M <- min_grid_pts + floor(sqrt(N))
  mseq <- seq_len(M)
  sM <- 1 - sqrt(M / (mseq - 0.5))
  Nflr <- floor(N / 4 + 0.5)
  b <- 1 / x[N] + sM / prior / x[Nflr]
  l <- N * lx(b, x)
  w <- 1 / vapply(mseq, FUN = function(j) sum(exp(l - l[j])), FUN.VALUE = 0)
  bdotw <- sum(b * w)
  k <- mean.default(log1p(-bdotw * x))
  sigma <- -k / bdotw
  if (wip) {
    k <- adjust_k_wip(k, n = N)
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
