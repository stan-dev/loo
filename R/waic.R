#' Widely applicable information criterion (WAIC)
#'
#' @export
#' @param log_lik an \eqn{S} by \eqn{N} matrix, where \eqn{S} is the size of the
#'   posterior sample (the number of simulations) and \eqn{N} is the number of
#'   data points. Typically (but not restricted to be) the object returned by
#'   \code{\link{extract_log_lik}}.
#'
#' @return a named list (of class \code{'loo'}) with components:
#'
#' \describe{
#' \item{\code{elpd_waic}}{expected log pointwise predictive density}
#' \item{\code{p_waic}}{estimated effective number of parameters}
#' \item{\code{waic}}{\code{-2 * elpd_waic} (i.e., converted to the deviance
#' scale)}
#' \item{\code{se_elpd_waic, se_p_waic, se_waic}}{standard errors}
#' \item{\code{pointwise}}{the pointwise contributions of each of the above
#' measures}
#' }
#'
#' @seealso \code{\link{loo-package}}, \code{\link{print.loo}},
#'   \code{\link{compare}}
#'
#' @examples
#' \dontrun{
#' log_lik1 <- extract_log_lik(stanfit1)
#' waic1 <- waic(log_lik1)
#' waic1
#'
#' log_lik2 <- extract_log_lik(stanfit2)
#' waic2 <- waic(log_lik2)
#' print(waic2, digits = 4)
#'
#' waic_diff <- compare(waic1, waic2)
#' print(waic_diff, digits = 2)
#' }
#'
waic <- function(log_lik) {
  if (!is.matrix(log_lik))
    stop('log_lik should be a matrix')
  pointwise <- pointwise_waic(log_lik)
  out <- totals(pointwise)
  nms <- names(pointwise)
  names(out) <- c(nms, paste0("se_", nms))
  out$pointwise <- cbind_list(pointwise)
  attr(out, "log_lik_dim") <- dim(log_lik)
  class(out) <- "loo"
  out
}
