.ic <- function(elpd) {
  stopifnot(is.vector(elpd), is.numeric(elpd))
  -2 * elpd
}

#' @importFrom matrixStats colVars
.pointwise_waic <- function(log_lik) {
  lpd <- logColMeansExp(log_lik)
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- .ic(elpd_waic)
  nlist(elpd_waic, p_waic, waic)
}
.pointwise_loo <- function(log_lik, vgis) {
  # vgis is output from vgisloo()
  lpd <- logColMeansExp(log_lik)
  elpd_loo <- vgis$loos
  p_loo <- lpd - elpd_loo
  looic <- .ic(elpd_loo)
  nlist(elpd_loo, p_loo, looic)
}
.totals <- function(pointwise) {
  N <- length(pointwise[[1L]])
  total  <- unlist_lapply(pointwise, sum)
  se <- sqrt(N * unlist_lapply(pointwise, var))
  as.list(c(total, se))
}
