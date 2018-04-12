#' Compute weighted expectations
#'
#' The \code{E_loo} function computes weighted expectations (means, variances,
#' quantiles) using the importance weights obtained from the
#' \link[=psis]{PSIS} smoothing procedure. The expectations estimated by the
#' \code{E_loo} function assume that the PSIS approximation is working well.
#' \strong{A small \link[=pareto-k-diagnostic]{Pareto k} estimate is necessary,
#' but not sufficient, for \code{E_loo} to give reliable estimates.} Additional
#' diagnostic checks for gauging the reliability of the estimates are in
#' development and will be added in a future release.
#'
#' @export
#' @param x A numeric vector or matrix.
#' @param psis_object An object returned by \code{\link{psis}}.
#' @param log_ratios Optionally, a vector or matrix (the same dimensions as
#'   \code{x}) of raw (not smoothed) log ratios. If working with log-likelihood
#'   values, the log ratios are the \strong{negative} of those values. If
#'   \code{log_ratios} is specified we are able to compute
#'   \link[=pareto-k-diagnostic]{Pareto k} diagnostics specific to \code{E_loo}.
#' @param type The type of expectation to compute. The options are
#'   \code{"mean"}, \code{"variance"}, and \code{"quantile"}.
#' @param probs For computing quantiles, a vector of probabilities.
#' @param ... Arguments passed to individual methods.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\code{value}}{
#'   The result of the computation.
#'
#'   For the matrix method, \code{value} is a vector with \code{ncol(x)}
#'   elements, with one exception: when \code{type} is \code{"quantile"} and
#'   multiple values are specified in \code{probs} the \code{value} component of
#'   the returned object is a \code{length(probs)} by \code{ncol(x)} matrix.
#'
#'   For the default/vector method the \code{value} component is scalar, with
#'   one exception: when \code{type} is \code{"quantile"} and multiple values
#'   are specified in \code{probs} the \code{value} component is a vector with
#'   \code{length(probs)} elements.
#'   }
#'  \item{\code{pareto_k}}{
#'   Function-specific diagnostic.
#'
#'   If \code{log_ratios} is not specified when calling \code{E_loo},
#'   \code{pareto_k} will be \code{NULL}. Otherwise, for the matrix method it
#'   will be a vector of length \code{ncol(x)} containing estimates of the shape
#'   parameter \eqn{k} of the generalized Pareto distribution. For the
#'   default/vector method, the estimate is a scalar.
#'  }
#' }
#'
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
#' log_ratios <- -1 * log_lik(fit)
#' dim(log_ratios)
#'
#' r_eff <- relative_eff(exp(-log_ratios), chain_id = rep(1:4, each = 1000))
#' psis_object <- psis(log_ratios, r_eff = r_eff, cores = 2)
#'
#' E_loo(yrep, psis_object, type = "mean")
#' E_loo(yrep, psis_object, type = "var")
#' E_loo(yrep, psis_object, type = "quantile", probs = 0.5) # median
#' E_loo(yrep, psis_object, type = "quantile", probs = c(0.1, 0.9))
#'
#' # To get Pareto k diagnostic with E_loo we also need to provide the negative
#' # log-likelihood values using the log_ratios argument.
#' E_loo(yrep, psis_object, type = "mean", log_ratios = log_ratios)
#' }
#'
E_loo <- function(x, psis_object, ...) {
  UseMethod("E_loo")
}

#' @rdname E_loo
#' @export
E_loo.default <-
  function(x,
           psis_object,
           ...,
           type = c("mean", "variance", "quantile"),
           probs = NULL,
           log_ratios = NULL) {
    stopifnot(
      is.numeric(x),
      is.psis(psis_object),
      length(x) == dim(psis_object)[1],
      is.null(log_ratios) || (length(x) == length(log_ratios))
    )
    type <- match.arg(type)
    E_fun <- .E_fun(type)
    r_eff <- NULL
    if (type == "variance") {
      r_eff <- attr(psis_object, "r_eff")
    }

    w <- as.vector(weights(psis_object, log = FALSE))
    x <- as.vector(x)
    out <- E_fun(x, w, probs, r_eff)

    if (is.null(log_ratios)) {
      warning("'log_ratios' not specified. Can't compute k-hat diagnostic.",
              call. = FALSE)
      khat <- NULL
    } else {
      khat <- E_loo_khat.default(x, psis_object, log_ratios)
    }
    list(value = out, pareto_k = khat)
  }

#' @rdname E_loo
#' @export
E_loo.matrix <-
  function(x,
           psis_object,
           ...,
           type = c("mean", "variance", "quantile"),
           probs = NULL,
           log_ratios = NULL) {
    stopifnot(
      is.numeric(x),
      is.psis(psis_object),
      identical(dim(x), dim(psis_object)),
      is.null(log_ratios) || identical(dim(x), dim(log_ratios))
    )
    type <- match.arg(type)
    E_fun <- .E_fun(type)
    fun_val <- numeric(1)
    r_eff <- NULL
    if (type == "variance") {
      r_eff <- attr(psis_object, "r_eff")
    } else if (type == "quantile") {
      stopifnot(
        is.numeric(probs),
        length(probs) >= 1,
        all(probs > 0 & probs < 1)
      )
      fun_val <- numeric(length(probs))
    }
    w <- weights(psis_object, log = FALSE)

    out <- vapply(seq_len(ncol(x)), function(i) {
      E_fun(x[, i], w[, i], probs = probs, r_eff = r_eff[i])
    }, FUN.VALUE = fun_val)

    if (is.null(log_ratios)) {
      warning("'log_ratios' not specified. Can't compute k-hat diagnostic.",
              call. = FALSE)
      khat <- NULL
    } else {
      khat <- E_loo_khat.matrix(x, psis_object, log_ratios)
    }
    list(value = out, pareto_k = khat)
  }



#' Select the function to use based on user's 'type' argument
#'
#' @noRd
#' @param type User's 'type' argument.
#' @return The function for computing the weighted expectation specified by
#'   'type'.
#'
.E_fun <- function(type = c("mean", "variance", "quantile")) {
  switch(
    type,
    "mean" = .wmean,
    "variance" = .wvar,
    "quantile" = .wquant
  )
}

#' loo-weighted mean, variance, and quantiles
#'
#' @noRd
#' @param x,w Vectors of the same length. This should be checked inside
#'   E_loo() before calling these functions.
#' @param probs Vector of probabilities.
#' @param ... ignored. Having ... allows 'probs' to be passed to .wmean and
#'   .wvar in E_loo() without resulting in an error.
#'
.wmean <- function(x, w, ...) {
  sum(w * x)
}
.wvar <- function(x, w, r_eff = NULL, ...) {
  if (is.null(r_eff)) {
    r_eff <- 1
  }
  r <- (x - .wmean(x, w))^2
  sum(w^2 * r) / r_eff
}
.wquant <- function(x, w, probs, ...) {
  if (all(w == w[1])) {
    return(quantile(x, probs = probs, names = FALSE))
  }

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

#' Compute function-specific k-hat diagnostics
#'
#' @noRd
#' @param log_ratios Vector or matrix of raw (not smoothed) log ratios with the
#'   same dimensions as \code{x}. If working with log-likelihood values, the log
#'   ratios are the negative of those values.
#' @return Vector (of length \code{NCOL(x)}) of k-hat estimates.
#'
E_loo_khat <- function(x, psis_object, log_ratios, ...) {
  UseMethod("E_loo_khat")
}

E_loo_khat.default <- function(x, psis_object, log_ratios, ...) {
  .E_loo_khat_i(x, log_ratios, attr(psis_object, "tail_len"))
}

E_loo_khat.matrix <- function(x, psis_object, log_ratios, ...) {
  tail_lengths <- attr(psis_object, "tail_len")
  sapply(seq_len(ncol(x)), function(i) {
    .E_loo_khat_i(x[, i], log_ratios[, i], tail_lengths[i])
  })
}

#' Compute function-specific khat estimates
#'
#' @noRd
#' @param x_i Vector of values of function h(theta)
#' @param log_ratios_i S-vector of log_ratios, log(r(theta)), for a single
#'   observation.
#' @param tail_len_i Integer tail length used for fitting GPD.
#' @return Scalar h-specific k-hat estimate.
#'
.E_loo_khat_i <- function(x_i, log_ratios_i, tail_len_i) {
    h_theta <- x_i
    r_theta <- exp(log_ratios_i - max(log_ratios_i))
    a <- sqrt(1 + h_theta^2) * r_theta
    log_a <- sort(log(a))

    S <- length(log_a)
    tail_ids <- seq(S - tail_len_i + 1, S)
    tail_sample <- log_a[tail_ids]
    cutoff <- log_a[min(tail_ids) - 1]

    smoothed <- psis_smooth_tail(tail_sample, cutoff)
    return(smoothed$k)
  }
