#' Widely applicable information criterion (WAIC)
#'
#' @export waic waic.array waic.matrix waic.function
#' @inheritParams loo
#' @param ... Currently ignored.
#'
#' @return A named list (of class \code{c("waic", "loo")}) with components:
#'
#' \describe{
#'  \item{\code{estimates}}{
#'  A matrix with two columns (\code{"Estimate"}, \code{"SE"}) and three
#'  rows (\code{"elpd_waic"}, \code{"p_waic"}, \code{"waic"}). This contains
#'  point estimates and standard errors of the expected log pointwise predictive
#'  density (\code{elpd_waic}), the effective number of parameters
#'  (\code{p_waic}) and the LOO information criterion \code{waic} (which is just
#'  \code{-2 * elpd_waic}, i.e., converted to deviance scale).
#'  }
#'  \item{\code{pointwise}}{
#'  A matrix with three columns (and number of rows equal to the number of
#'  observations) containing the pointwise contributions of each of the above
#'  measures (\code{elpd_waic}, \code{p_waic}, \code{waic}).
#'  }
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
  pointwise <- pointwise_waic(log_lik = ll)
  throw_pwaic_warnings(out$pointwise[, "p_waic"], digits = 1)
  waic_object(pointwise, lldim)
}

#' @export
#' @templateVar fn waic
#' @template function
#'
waic.function <- function(x, args, ...) {
  lldim <- with(args, c(S, N))
  pointwise <- pointwise_waic(llfun = x, llargs = args)
  throw_pwaic_warnings(out$pointwise[, "p_waic"], digits = 1)
  waic_object(pointwise, lldim)
}



# internal ----------------------------------------------------------------

# structure the object returned by the waic methods
waic_object <- function(pointwise, log_lik_dim) {
  estimates <- table_of_estimates(pointwise)
  structure(
    nlist(estimates, pointwise),
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


pointwise_waic <- function(log_lik, llfun = NULL, llargs = NULL) {
  if (!missing(log_lik)) {
    lpd <- logColMeansExp(log_lik)
    p_waic <- matrixStats::colVars(log_lik)
  } else {
    if (is.null(llfun) || is.null(llargs))
      stop("Either 'log_lik' or 'llfun' and 'llargs' must be specified.",
           call. = FALSE)
    lpd <- logColMeansExp_llfun(llfun, llargs)
    p_waic <- colVars_llfun(llfun, llargs)
  }
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  cbind(elpd_waic, p_waic, waic)
}

colVars_llfun <- function(fun, args) {
  vapply(seq_len(args$N), FUN = function(i) {
    x <- fun(i = i, data = args$data[i,,drop=FALSE], draws = args$draws)
    var(as.vector(x))
  }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
}

