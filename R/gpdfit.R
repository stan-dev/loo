#' Generalized Pareto distribution
#'
#' Estimate the parameters \eqn{k} and \eqn{\sigma} of the generalized Pareto
#' distribution, given a sample \eqn{x}.
#'
#' @keywords internal
#' @export
#' @param x A numeric vector. The sample from which to estimate the parameters.
#' @return A named list with components \code{k} and \code{sigma}.
#'
#' @details Here the parameter \eqn{k} is the negative of \eqn{k} in Zhang &
#'   Stephens (2009).
#'
#' @note This function is primarily intended for internal use, but is exported
#'   so that users can call it directly if desired. Users simply wishing to
#'   compute LOO should use the \code{\link{loo}} function.
#'
#' @seealso \code{\link{psislw}}, \code{\link{loo}},
#' \code{\link{loo-package}}
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
  M <- 80 + floor(sqrt(N))  # note: original paper used  20 + floor(sqrt(n))
  mseq <- seq_len(M)
  sM <- 1 - sqrt(M / (mseq - 0.5))
  Nflr <- floor(N / 4 + 0.5)
  b <- 1 / x[N] + sM / prior / x[Nflr]
  l <- N * lx(b, x)
  w <- 1 / vapply(mseq, FUN = function(j) sum(exp(l - l[j])), FUN.VALUE = 0)
  bdotw <- sum(b * w)
  k <- mean.default(log1p(-bdotw * x))
  sigma <- -k / bdotw
  nlist(k, sigma)
}
