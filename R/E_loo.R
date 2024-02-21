#' Compute weighted expectations
#'
#' The `E_loo()` function computes weighted expectations (means, variances,
#' quantiles) using the importance weights obtained from the
#' [PSIS][psis()] smoothing procedure. The expectations estimated by the
#' `E_loo()` function assume that the PSIS approximation is working well.
#' **A small [Pareto k][pareto-k-diagnostic] estimate is necessary,
#' but not sufficient, for `E_loo()` to give reliable estimates.** Additional
#' diagnostic checks for gauging the reliability of the estimates are in
#' development and will be added in a future release.
#'
#' @export
#' @param x A numeric vector or matrix.
#' @param psis_object An object returned by [psis()].
#' @param log_ratios Optionally, a vector or matrix (the same dimensions as `x`)
#'   of raw (not smoothed) log ratios. If working with log-likelihood values,
#'   the log ratios are the **negative** of those values. If `log_ratios` is
#'   specified we are able to compute more accurate [Pareto k][pareto-k-diagnostic]
#'   diagnostics specific to `E_loo()`.
#' @param type The type of expectation to compute. The options are
#'   `"mean"`, `"variance"`, `"sd"`, and `"quantile"`.
#' @param probs For computing quantiles, a vector of probabilities.
#' @param ... Arguments passed to individual methods.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`value`}{
#'   The result of the computation.
#'
#'   For the matrix method, `value` is a vector with `ncol(x)`
#'   elements, with one exception: when `type="quantile"` and
#'   multiple values are specified in `probs` the `value` component of
#'   the returned object is a `length(probs)` by `ncol(x)` matrix.
#'
#'   For the default/vector method the `value` component is scalar, with
#'   one exception: when `type="quantile"` and multiple values
#'   are specified in `probs` the `value` component is a vector with
#'   `length(probs)` elements.
#'   }
#'  \item{`pareto_k`}{
#'   Function-specific diagnostic.
#'
#'   For the matrix method it will be a vector of length `ncol(x)`
#'   containing estimates of the shape parameter \eqn{k} of the
#'   generalized Pareto distribution. For the default/vector method,
#'   the estimate is a scalar. If `log_ratios` is not specified when
#'   calling `E_loo()`, the smoothed log-weights are used to estimate
#'   Pareto-k's, which may produce optimistic estimates.
#'
#'   For `type="mean"`, `type="var"`, and `type="sd"`, the returned Pareto-k is
#'   the maximum of the Pareto-k's for the left and right tail of \eqn{hr} and
#'   the right tail of \eqn{r}, where \eqn{r} is the importance ratio and
#'   \eqn{h=x} for `type="mean"` and \eqn{h=x^2} for `type="var"` and
#'   `type="sd"`. For `type="quantile"`, the returned Pareto-k is the Pareto-k
#'   for the right tail of \eqn{r}.
#'  }
#' }
#'
#'
#' @examples
#' \donttest{
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
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
#' fit <- stan_glm(weight ~ group, data = d, refresh = 0)
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
#' E_loo(yrep, psis_object, type = "sd")
#' E_loo(yrep, psis_object, type = "quantile", probs = 0.5) # median
#' E_loo(yrep, psis_object, type = "quantile", probs = c(0.1, 0.9))
#'
#' # We can get more accurate Pareto k diagnostic if we also provide
#' # the log_ratios argument
#' E_loo(yrep, psis_object, type = "mean", log_ratios = log_ratios)
#' }
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
           type = c("mean", "variance", "sd", "quantile"),
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
    w <- as.vector(weights(psis_object, log = FALSE))
    x <- as.vector(x)
    out <- E_fun(x, w, probs)

    if (is.null(log_ratios)) {
      # Use of smoothed ratios gives slightly optimistic
      # Pareto-k's, but these are still better than nothing
      log_ratios <- weights(psis_object, log = TRUE)
    }
    h <- switch(
      type,
      "mean" = x,
      "variance" = x^2,
      "sd" = x^2,
      "quantile" = NULL
    )
    khat <- E_loo_khat.default(h, psis_object, log_ratios)
    list(value = out, pareto_k = khat)
  }

#' @rdname E_loo
#' @export
E_loo.matrix <-
  function(x,
           psis_object,
           ...,
           type = c("mean", "variance", "sd", "quantile"),
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
    if (type == "quantile") {
      stopifnot(
        is.numeric(probs),
        length(probs) >= 1,
        all(probs > 0 & probs < 1)
      )
      fun_val <- numeric(length(probs))
    }
    w <- weights(psis_object, log = FALSE)

    out <- vapply(seq_len(ncol(x)), function(i) {
      E_fun(x[, i], w[, i], probs = probs)
    }, FUN.VALUE = fun_val)

    if (is.null(log_ratios)) {
      # Use of smoothed ratios gives slightly optimistic
      # Pareto-k's, but these are still better than nothing
      log_ratios <- weights(psis_object, log = TRUE)
    }
    h <- switch(
      type,
      "mean" = x,
      "variance" = x^2,
      "sd" = x^2,
      "quantile" = NULL
    )
    khat <- E_loo_khat.matrix(h, psis_object, log_ratios)
    list(value = out, pareto_k = khat)
  }



#' Select the function to use based on user's 'type' argument
#'
#' @noRd
#' @param type User's `type` argument.
#' @return The function for computing the weighted expectation specified by
#'   `type`.
#'
.E_fun <- function(type = c("mean", "variance", "sd", "quantile")) {
  switch(
    type,
    "mean" = .wmean,
    "variance" = .wvar,
    "sd" = .wsd,
    "quantile" = .wquant
  )
}

#' loo-weighted mean, variance, and quantiles
#'
#' @noRd
#' @param x,w Vectors of the same length. This should be checked inside
#'   `E_loo()` before calling these functions.
#' @param probs Vector of probabilities.
#' @param ... ignored. Having ... allows `probs` to be passed to `.wmean()` and
#'   `.wvar()` in `E_loo()` without resulting in an error.
#'
.wmean <- function(x, w, ...) {
  sum(w * x)
}
.wvar <- function(x, w, ...) {
  # The denominator (1- sum(w^2)) is equal to (ESS-1)/ESS, where effective
  # sample size ESS is estimated with the generic target quantity invariant
  # estimate 1/sum(w^2), see e.g. "Monte Carlo theory, methods and examples"
  # by Owen (2013).
  (sum(.wmean(x^2, w)) - sum(.wmean(x, w)^2)) / (1 - sum(w^2))
}
.wsd <- function(x, w, ...) {
  sqrt(.wvar(x, w))
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
#'   same dimensions as `x`. If working with log-likelihood values, the log
#'   ratios are the negative of those values.
#' @return Vector (of length `NCOL(x)`) of k-hat estimates.
#'
E_loo_khat <- function(x, psis_object, log_ratios, ...) {
  UseMethod("E_loo_khat")
}
#' @export
E_loo_khat.default <- function(x, psis_object, log_ratios, ...) {
  .E_loo_khat_i(x, log_ratios, attr(psis_object, "tail_len"))
}
#' @export
E_loo_khat.matrix <- function(x, psis_object, log_ratios, ...) {
  tail_lengths <- attr(psis_object, "tail_len")
  if (is.null(x)) {
    sapply(seq_len(ncol(log_ratios)), function(i) {
      .E_loo_khat_i(x, log_ratios[, i], tail_lengths[i])
    })
  } else {
    sapply(seq_len(ncol(log_ratios)), function(i) {
      .E_loo_khat_i(x[, i], log_ratios[, i], tail_lengths[i])
    })
  }
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
  khat_r <- posterior::pareto_khat(r_theta, tail = "right", ndraws_tail = tail_len_i)$khat
  if (is.null(x_i)) {
    khat_r
  } else {
    khat_hr <- posterior::pareto_khat(h_theta * r_theta, tail = "both", ndraws_tail = tail_len_i)$khat
    max(khat_hr, khat_r)
  }
}
