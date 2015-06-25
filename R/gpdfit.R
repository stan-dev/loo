#' Generalized Pareto distribution
#'
#' Estimate the parameters \eqn{k} and \eqn{\sigma} of the generalized Pareto
#' distribution, given a sample \eqn{x}.
#'
#' @keywords internal
#' @export
#' @param x a numeric vector. The sample from which to estimate the parameters.
#' @return a list. Estimates for \eqn{k} and \eqn{\sigma}.
#'
#' @details Here the parameter \eqn{k} is the negative of \eqn{k} in Zhang &
#'   Stephens (2009).
#'
#' @note This function is primarily intended for internal use, but is exported
#'   so that users can call it directly if desired. Users simply wishing to
#'   compute LOO and WAIC should use the \code{\link{loo_and_waic}} function.
#'
#' @seealso \code{\link{vgislw}}, \code{\link{loo_and_waic}},
#' \code{\link{loo-package}}
#'
#' @references
#' Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
#' for the generalized Pareto distribution. \emph{Technometrics} \strong{51},
#' 316-325.
#'

gpdfit <- function(x) {
  n <- length(x)
  x <- sort.int(x, method = "quick")
  prior <- 3
  m <- 80 + floor(sqrt(n))  # note: original paper used  20 + floor(sqrt(n))
  b <- 1 / x[n] + (1 - sqrt(m / seq_min_half(m))) / prior / x[floor(n / 4 + 0.5)]
  L <- vapply_seq(m, function(i) n * lx(b[i], x))
  w <- 1 / vapply_seq(m, function(i) sum(exp(L - L[i])))
  b <- sum(b * w)
  k <- mean.default(log(1 - b * x))
  sigma <- -k / b
  nlist(k, sigma)
}
