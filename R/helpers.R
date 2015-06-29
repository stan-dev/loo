# waic and loo helpers ----------------------------------------------------
#' @importFrom matrixStats colLogSumExps
logColMeansExp <- function(x) {
  # should be more stable than log(colMeans(exp(x)))
  S <- nrow(x)
  colLogSumExps(x) - log(S)
}

#' @importFrom matrixStats colVars
pointwise_waic <- function(log_lik) {
  lpd <- logColMeansExp(log_lik)
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  nlist(elpd_waic, p_waic, waic)
}
pointwise_loo <- function(log_lik, vgis) {
  # vgis is output from vgisloo()
  lpd <- logColMeansExp(log_lik)
  elpd_loo <- vgis$loos
  p_loo <- lpd - elpd_loo
  looic <- -2 * elpd_loo
  nlist(elpd_loo, p_loo, looic)
}
totals <- function(pointwise) {
  N <- length(pointwise[[1L]])
  total  <- unlist_lapply(pointwise, sum)
  se <- sqrt(N * unlist_lapply(pointwise, var))
  as.list(c(total, se))
}


# VGIS helpers ------------------------------------------------------------

# inverse-CDF of generalized Pareto distribution (formula from Wikipedia)
qgpd <- function(p, xi = 1, mu = 0, beta = 1, lower.tail = TRUE) {
  if (!lower.tail)
    p <- 1 - p
  mu + beta * ((1 - p)^(-xi) - 1) / xi
}

# lx <- function(a, x) {
#   k <- mean.default(log1p(-a * x))
#   log(-a / k) - k - 1
# }
lx <- function(a, x) {
  # vectorized version of old lx
  b <- -a
  bx <- outer(x, b)
  d <- dim(bx)
  k <- .colMeans(log1p(bx), d[1], d[2])
  log(b / k) - k - 1
}

#' @importFrom matrixStats logSumExp
lw_normalize <- function(y) {
  y - logSumExp(y)
}
lw_truncate <- function(y, wtrunc) {
  if (wtrunc == 0)
    return(y)
  logS <- log(length(y))
  lwtrunc <- wtrunc * logS - logS + logSumExp(y)
  y[y > lwtrunc] <- lwtrunc
  y
}
lw_cutpoint <- function(y, wcp, min_cut) {
  cp <- quantile(y, 1 - wcp, names = FALSE)
  max(cp, min_cut)
}

# The parallelization functions mclapply and parLapply return a list of lists:
# vgis is a list of length N=ncol(lw). Each of the N elements of vgis is itself
# a list of length 2. In each of these N lists of length 2 the first component
# is a vector of length S=nrow(lw) containing the modified log weights and the
# second component is the estimate of the pareto shape parameter k. ux is now a
# list of length 2*N. the odd elements contain the modified log weights and the
# even elements contain the pareto k estimates. This function cbinds the log
# weight vectors into a matrix and concatenates the k estimates into a vector.
.vgis_out <- function(vgis) {
  L <- length(vgis)
  ux <- unlist(vgis, recursive = FALSE, use.names = FALSE)
  lw_smooth <- cbind_list(ux[nodds(L)])
  pareto_k <- unlist(ux[nevens(L)])
  nlist(lw_smooth, pareto_k)
}
