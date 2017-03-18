#' Compute weighted expectations
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
loo_expectation <- function(x, lw, ...) {
  UseMethod("loo_expectation")
}

#' @rdname loo_expectation
#' @export
loo_expectation.default <-
  function(x,
           lw,
           ...,
           type = c("mean", "var", "quantile"),
           probs) {
    stopifnot(is.numeric(x), is.numeric(lw), length(lw) == length(x))
    E_loo <- .E_loo(type)
    x <- as.vector(x)
    w <- exp(as.vector(lw))
    E_loo(x, w, probs)
  }

#' @rdname loo_expectation
#' @export
loo_expectation.matrix <-
  function(x,
           lw,
           ...,
           type = c("mean", "var", "quantile"),
           probs) {
    stopifnot(is.numeric(x), is.numeric(lw), identical(dim(x), dim(lw)))
    type <- match.arg(type)
    E_loo <- .E_loo(type)
    fun_val <- if (type == "quantile")
      numeric(length(probs)) else numeric(1)

    w <- exp(lw)
    vapply(seq_len(ncol(x)), function(k) {
      E_loo(x[, k], w[, k], probs)
    }, FUN.VALUE = fun_val)
  }



# @param type user's 'type' argument
# @return the function for computing the expectation specified by 'type'
.E_loo <- function(type = c("mean", "var", "quantile")) {
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
#   loo_expectation() before calling these functions.
# @param probs vector of probabilities.
# @param ... ignored. having ... allows 'probs' to be passed to .wmean and .wvar
#   in loo_expectation() without resulting in an error.
#
.wmean <- function(x, w, ...) {
  sum(w * x)
}
.wvar <- function(x, w, ...) {
  r <- (x - .wmean(x, w))^2
  sum(w * r)
}
.wquant <- function(x, w, probs) {
  stopifnot(all(probs > 0 | probs < 1))
  x <- sort(x)
  ww <- cumsum(w)
  ww <- ww / ww[length(ww)]

  y <- numeric(length(probs))
  for (j in seq_along(probs)) {
    ids <- which(ww >= probs[j])
    wi <- min(ids)
    if (wi == 1) {
      y[j] <- x[1]
    } else {
      w1 <- ww[wi - 1]
      x1 <- x[wi - 1]
      y[j] <- x1 + (x[wi] - x1) * (probs[j] - w1) / (ww[wi] - w1)
    }
  }
  return(y)
}
