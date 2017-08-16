#' Generalized Pareto distribution
#'
#' Estimate the parameters \eqn{k} and \eqn{\sigma} of the generalized Pareto
#' distribution, given a sample \eqn{x}. The Pareto fit uses a prior for
#' \eqn{k}, which will stabilize estimates for very small Monte Carlo sample
#' sizes and low effective sample sizes. The weakly informative prior is a
#' Gaussian prior centered on 0.5.
#'
#' @keywords internal
#' @export
#' @param x A numeric vector. The sample from which to estimate the parameters.
#' @return A named list with components \code{k} and \code{sigma}.
#'
#' @details Here the parameter \eqn{k} is the negative of \eqn{k} in Zhang &
#'   Stephens (2009).
#'
#' @template internal-function-note
#'
#' @seealso \code{\link{psis}}, \code{\link{pareto-k-diagnostic}},
#'   \code{\link{loo-package}}
#'
#' @references
#' Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
#' for the generalized Pareto distribution. \emph{Technometrics} \strong{51},
#' 316-325.
#'
gpdfit <- function(x) {
  N <- length(x)
  x <- sort.int(x, method = "quick")
  prior <- 3
  M <- 80 + floor(sqrt(N))
  mseq <- seq_len(M)
  sM <- 1 - sqrt(M / (mseq - 0.5))
  Nflr <- floor(N / 4 + 0.5)
  b <- 1 / x[N] + sM / prior / x[Nflr]
  l <- N * lx(b, x)
  w <- 1 / vapply(mseq, FUN = function(j) sum(exp(l - l[j])), FUN.VALUE = 0)
  bdotw <- sum(b * w)
  k <- mean.default(log1p(-bdotw * x))
  sigma <- -k / bdotw
  k <- adjust_k(k, n = N)
  nlist(k, sigma)
}


# internal ----------------------------------------------------------------

lx <- function(a,x) {
  a <- -a
  k <- sapply(a, FUN = function(y) mean(log1p(y * x)))
  log(a / k) - k - 1
}

# Adjust k based on weakly informative Gaussian prior centered on 0.5. This will
# stabilize estimates for very small Monte Carlo sample sizes and low neff
# cases.
# @param k khat estimate
# @param n number of tail samples used to fit GPD
adjust_k <- function(k, n) {
  a <- 10
  nplusa <- n + a
  k * n / nplusa + a * 0.5 / nplusa
}
