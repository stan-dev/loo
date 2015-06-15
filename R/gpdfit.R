#' Estimate the parameters \eqn{k} and \eqn{\sigma} of the generalized Pareto
#' distribution, given a sample \eqn{x}.
#'
#' @export
#' @param x a numeric vector. The sample from which to estimate the parameters.
#' @return a list. Estimates for \eqn{k} and \eqn{\sigma}.
#' @note Here the parameter \eqn{k} is the negative of the \eqn{k} in the paper
#' of Zhang & Stephens.
#' @references Zhang & Stephens (2009)
#'

gpdfit <- function(x) {
  n <- length(x)
  x <- sort.int(x)
  prior <- 3
  m <- 80 + floor(sqrt(n))  # note: original paper used m <- 20+floor(sqrt(n))
  b <- 1/x[n] + (1 - sqrt(m/seq_min_half(m)))/prior/x[floor(n/4 + 0.5)]
  L <- vapply1m(m, function(i) n * lx(b[i], x))
  w <- 1 / vapply1m(m, function(i) sum(exp(L - L[i])))
  b <- sum(b * w)
  k <- mean.default(log(1 - b * x))
  sigma <- -k/b
  nlist(k, sigma)
}
