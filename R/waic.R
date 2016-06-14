#' Widely applicable information criterion (WAIC)
#'
#' @export waic waic.matrix waic.function
#' @inheritParams loo
#' @param ... Other arguments. Currently ignored.
#'
#' @return A named list (of class \code{'loo'}) with components:
#'
#' \describe{
#' \item{\code{elpd_waic, se_elpd_waic}}{expected log pointwise predictive
#' density and standard error}
#' \item{\code{p_waic, se_p_waic}}{estimated effective number of parameters and
#' standard error}
#' \item{\code{waic, se_waic}}{\code{-2 * elpd_waic} (i.e., converted to the
#' deviance scale) and standard error}
#' \item{\code{pointwise}}{the pointwise contributions of each of the above
#' measures}
#' }
#'
#' @seealso \code{\link{compare}}, \code{\link{print.loo}},
#' \code{\link{loo-package}}
#'
#'
#' @examples
#' \dontrun{
#' log_lik1 <- extract_log_lik(stanfit1)
#' log_lik2 <- extract_log_lik(stanfit2)
#' (waic1 <- waic(log_lik1))
#' (waic2 <- waic(log_lik2))
#' print(compare(waic1, waic2), digits = 2)
#' }
#'
waic <- function(x, ...) {
  UseMethod("waic")
}

#' @export
#' @templateVar fn waic
#' @template matrix
#'
waic.matrix <- function(x, ...) {
  if (any(is.na(x)))
    stop("NA log-likelihood values found.")
  out <- pointwise_waic(log_lik = x)
  pwaic_warnings(out$pointwise[, "p_waic"], digits = 1)
  structure(out, log_lik_dim = dim(x), class = "loo")
}

#' @export
#' @templateVar fn waic
#' @template function
#'
waic.function <- function(x, ..., args) {
  if (missing(args))
    stop("'args' must be specified.")
  out <- pointwise_waic(llfun = x, llargs = args)
  pwaic_warnings(out$pointwise[, "p_waic"], digits = 1)
  structure(out, log_lik_dim = with(args, c(S,N)), class = "loo")
}
