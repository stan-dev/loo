#' Compute weighted expectations
#'
#' The \code{E_loo} function computes weighted expectations (means, variances,
#' quantiles) using the smoothed importance weights obtained from the
#' \link[=psislw]{PSIS} procedure. The expectations estimated by the
#' \code{E_loo} function assume that the PSIS approximation is working well.
#' \strong{A small \link[=pareto-k-diagnostic]{Pareto k} estimate is necessary,
#' but not sufficient, for \code{E_loo} to give reliable estimates.} Additional
#' diagnostic checks for gauging the reliability of the estimates are in
#' development and will be added in a future release.
#'
#' @export
#' @param x A numeric matrix or vector.
#' @param lw A numeric matrix (or vector) of smoothed log-weights with the same
#'   dimensions (or length) as \code{x}. Typically this will be the
#'   \code{"lw_smooth"} element of the list returned by \code{\link{psislw}}.
#' @param type The type of expectation to compute. The options are
#'   \code{"mean"}, \code{"var"} (variance), and \code{"quantile"}.
#' @param probs A vector of probabilities. Ignored unless \code{type} is
#'   \code{"quantile"}.
#' @param ... For the generic function, arguments passed to the
#'   individual methods.
#'
#' @return The matrix method returns a vector with \code{ncol(x)} elements, with
#'   one exception: when \code{type} is \code{"quantile"} and multiple values
#'   are specified in \code{probs} the returned object is a \code{length(probs)}
#'   by \code{ncol(x)} matrix.
#'
#'   The default/vector method returns a scalar, with one exception: when
#'   \code{type} is \code{"quantile"} and multiple values are specified in
#'   \code{probs} the returned object is a vector with \code{length(probs)}
#'   elements.
#'
#' @examples
#' \donttest{
#' # Use rstanarm package to quickly fit a model and get both a log-likelihood
#' # matrix and draws from the posterior predictive distribution
#' library("rstanarm")
#'
#' # data from help("lm")
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' d <- data.frame(
#'   weight = c(ctl, trt),
#'   group = gl(2, 10, 20, labels = c("Ctl","Trt"))
#' )
#' fit <- stan_glm(weight ~ group, data = d)
#' yrep <- posterior_predict(fit)
#' dim(yrep)
#'
#' ll <- log_lik(fit)
#' lw <- psislw(-ll, cores = 2)$lw_smooth
#' dim(lw)
#'
#' E_loo(yrep, lw, type = "mean")
#' E_loo(yrep, lw, type = "var")
#' E_loo(yrep, lw, type = "quantile", probs = 0.5) # median
#' E_loo(yrep, lw, type = "quantile", probs = c(0.1, 0.9))
#' }
#'
E_loo <- function(x, lw, ...) {
  UseMethod("E_loo")
}

#' @rdname E_loo
#' @export
E_loo.default <-
  function(x,
           lw,
           ...,
           type = c("mean", "var", "quantile"),
           probs) {
    stopifnot(is.numeric(x), is.numeric(lw), length(lw) == length(x))
    E_fun <- .E_fun(type)
    x <- as.vector(x)
    w <- exp(as.vector(lw))
    E_fun(x, w, probs)
  }

#' @rdname E_loo
#' @export
E_loo.matrix <-
  function(x,
           lw,
           ...,
           type = c("mean", "var", "quantile"),
           probs) {
    stopifnot(is.numeric(x), is.numeric(lw), identical(dim(x), dim(lw)))
    type <- match.arg(type)
    E_fun <- .E_fun(type)
    if (type == "quantile") {
      stopifnot(is.numeric(probs), length(probs) >= 1)
      fun_val <- numeric(length(probs))
    } else {
      fun_val <- numeric(1)
    }

    w <- exp(lw)
    vapply(seq_len(ncol(x)), function(i) {
      E_fun(x[, i], w[, i], probs)
    }, FUN.VALUE = fun_val)
  }



# @param type user's 'type' argument
# @return the function for computing the expectation specified by 'type'
.E_fun <- function(type = c("mean", "var", "quantile")) {
  switch(
    match.arg(type),
    "mean" = .wmean,
    "var" = .wvar,
    "quantile" = .wquant
  )
}

# loo-weighted mean, variance, and quantiles
#
# @param x,w vectors of the same length. this should be checked inside
#   E_loo() before calling these functions.
# @param probs vector of probabilities.
# @param ... ignored. having ... allows 'probs' to be passed to .wmean and .wvar
#   in E_loo() without resulting in an error.
#
.wmean <- function(x, w, ...) {
  sum(w * x)
}
.wvar <- function(x, w, ...) {
  r <- (x - .wmean(x, w))^2
  sum(w * r)
}
.wquant <- function(x, w, probs) {
  stopifnot(all(probs > 0 & probs < 1))
  if (all(w == w[1]))
    return(quantile(x, probs = probs, names = FALSE))

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  ww <- cumsum(w)
  ww <- ww / ww[length(ww)]

  qq <- numeric(length(probs))
  for (j in seq_along(probs)) {
    ids <- which(ww >= probs[j])
    wi <- min(ids)
    if (wi == 1) {
      qq[j] <- x[1]
    } else {
      w1 <- ww[wi - 1]
      x1 <- x[wi - 1]
      qq[j] <- x1 + (x[wi] - x1) * (probs[j] - w1) / (ww[wi] - w1)
    }
  }
  return(qq)
}
