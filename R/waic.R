#' Widely applicable information criterion (WAIC)
#'
#' @export waic waic.array waic.matrix waic.function
#' @inheritParams loo
#' @param ... Other arguments. Currently ignored.
#'
#' @return A named list (of class \code{c("waic", "loo")}) with components:
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
#' @template array
#'
waic.array <- function(x, ...) {
  waic.matrix(llarray_to_matrix(x), ...)
}

#' @export
#' @templateVar fn waic
#' @template matrix
#'
waic.matrix <- function(x, ...) {
  ll <- validate_ll(x)
  lldim <- dim(ll)
  out <- pointwise_waic(log_lik = ll)
  throw_pwaic_warnings(out$pointwise[, "p_waic"], digits = 1)
  waic_object(out, lldim)
}

#' @export
#' @templateVar fn waic
#' @template function
#'
waic.function <- function(x, args, ...) {
  lldim <- with(args, c(S, N))
  out <- pointwise_waic(llfun = x, llargs = args)
  throw_pwaic_warnings(out$pointwise[, "p_waic"], digits = 1)
  waic_object(out, lldim)
}



# internal ----------------------------------------------------------------

# structure the object returned by the waic methods
waic_object <- function(object, log_lik_dim) {
  structure(
    object,
    log_lik_dim = log_lik_dim,
    class = c("waic", "loo")
  )
}

# waic warnings
# @param p 'p_waic' estimates
throw_pwaic_warnings <- function(p, digits = 1) {
  badp <- p > 0.4
  if (any(badp)) {
    count <- sum(badp)
    prop <- count / length(badp)
    .warn(paste0(count, " (", .fr(100 * prop, digits),
                 "%) p_waic estimates greater than 0.4."),
          "\nWe recommend trying loo() instead.")
  }
  invisible(NULL)
}

