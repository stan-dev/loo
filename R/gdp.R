#' The Generalized Pareto Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the generalized Pareto distribution with location equal to \code{mu},
#' scale equal to \code{sigma}, and shape equal to \code{k}.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#' to be the number required.
#' @param mu scalar location parameter
#' @param sigma scalar, positive scale parameter
#' @param k scalar shape parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#' otherwise, \eqn{P[X > x]}.
#'
#' @name GPD


#' @rdname GPD
#' @export
dgpd <- function(x, mu = 0, sigma = 1, k = 0, log = FALSE) {
  stopifnot(length(mu) == 1 && length(sigma) == 1 && length(k) == 1)
  if (is.na(sigma) || sigma <= 0) {
    return(rep(NaN, length(x)))
  }
  d <- (x - mu) / sigma
  ind <- (d > 0) & ((1 + k * d) > 0)
  ind[is.na(ind)] <- FALSE
  if (k == 0) {
    d[ind] <- -d[ind] - log(sigma)
  } else {
    d[ind] <- log1p(k * d[ind]) * -(1 / k  + 1) - log(sigma)
  }
  d[!ind] <- -Inf
  if (!log) {
    d <- exp(d)
  }
  d
}


#' @rdname GPD
#' @export
pgpd <- function(q, mu = 0, sigma = 1, k = 0, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(length(mu) == 1 && length(sigma) == 1 && length(k) == 1)
  if (is.na(sigma) || sigma <= 0) {
    return(rep(NaN, length(q)))
  }
  q <- pmax(q - mu, 0) / sigma
  if (k == 0) {
    p <- 1 - exp(-q)
  } else {
    p <- -expm1(log(pmax(1 + k * q, 0)) * -(1 / k))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

#' @rdname GPD
#' @export
qgpd <- function(p, mu = 0, sigma = 1, k = 0, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(length(mu) == 1 && length(sigma) == 1 && length(k) == 1)
  if (is.na(sigma) || sigma <= 0) {
    return(rep(NaN, length(p)))
  }
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (k == 0) {
    q <-  mu - sigma * log1p(-p)
  } else {
    q <- mu + sigma * expm1(-k * log1p(-p)) / k
  }
  q
}

#' @rdname GPD
#' @export
rgpd <- function(n, mu = 0, sigma = 1, k = 0) {
  stopifnot(
    length(n) == 1 && length(mu) == 1 && length(sigma) == 1 && length(k) == 1
  )
  if (is.na(sigma) || sigma <= 0) {
    return(rep(NaN, n))
  }
  if (k == 0) {
    r <- mu + sigma * rexp(n)
  } else {
    r <- mu + sigma * expm1(-k * log(runif(n))) / k
  }
  r
}
