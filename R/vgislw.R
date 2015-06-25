#' VGIS: Very good importance sampling
#'
#' @export
#' @param lw a matrix or vector of log weights. For for computing LOO \code{lw =
#'   -log_lik} (see \code{\link{extract_log_lik}}) and is an \eqn{S} by \eqn{N}
#'   matrix where \eqn{S} is the number of simulations and \eqn{N} is the number
#'   of data points. (If \code{lw} is a vector it will be coerced to a
#'   one-column matrix.)
#' @param wcp the proportion of importance weights to use for the generalized
#'   Pareto fit. The \code{100*wcp}\% largest weights are used as the sample
#'   from which to estimate the parameters of the generalized Pareto
#'   distribution.
#' @param thresh when used for the generalized pareto fit, the largest
#'   \code{100*wcp}\% of the importance weights are modified according to
#'   \code{pmax(x, max(x) - thresh)} for numerical stability.
#' @param kmax maximum allowed value for the generalized Pareto shape parameter
#'   \eqn{k}.
#' @param wtrunc for truncating very large weights to \eqn{S}^\code{wtrunc}. Set
#'   to zero for no truncation.
#' @param cores the number of cores to use for parallelization.
#'
#' @return A named with list with components \code{lw_smooth} (modified log
#'   weights) and \code{pareto_k} (estimated generalized Pareto shape parameters
#'   \eqn{k}).
#'
#'
#' @details See the 'VGIS-LOO' section in \code{\link{loo-package}}.
#'
#' @note This function is primarily intended for internal use, but is exported
#'   so that users can call it directly for other purposes. Users simply
#'   wishing to compute LOO and WAIC should use the \code{\link{loo_and_waic}}
#'   function.
#'
#' @seealso \code{\link{loo-package}}, \code{\link{loo_and_waic}}
#'
#' @importFrom matrixStats logSumExp
#' @importFrom parallel mclapply makePSOCKcluster stopCluster parLapply
#'
vgislw <- function(lw, wcp = 0.2, thresh = 100, kmax = 2, wtrunc = 3/4,
                   cores = parallel::detectCores()) {
  .vgis <- function(i) {
    x <- lw[, i]
    S <- length(x)
    # split into body and right tail
    cutoff <- quantile(x, 1 - wcp, names = FALSE)
    x_cut <- x > cutoff
    xbody <- x[!x_cut]
    xtail <- x[x_cut]
    tail_len <- length(xtail)
    # store order of tail samples
    tail_ord <- order(xtail)
    # fit generalized Pareto distribution to the right tail samples
    xtail <- pmax(xtail, max(xtail) - thresh)
    exp_cutoff <- exp(cutoff)
    fit <- gpdfit(exp(xtail) - exp_cutoff)
    k <- min(fit$k, kmax)
    # compute order statistics for the fit
    qq <- qgpd(seq_min_half(tail_len)/tail_len, xi = k, beta = fit$sigma)
    qq <- qq + exp_cutoff
    # remap back to the original order
    slq <- rep.int(0, tail_len)
    slq[tail_ord] <- log(qq)
    # join body and gPd smoothed tail
    qx <- x
    qx[!x_cut] <- xbody
    qx[x_cut] <- slq
    if (wtrunc > 0) {
      # truncate
      logS <- log(S)
      lwtrunc <- wtrunc * logS - logS + logSumExp(qx)
      qx[qx > lwtrunc] <- lwtrunc
    }
    # renormalize weights
    lwx <- qx - logSumExp(qx)
    # return log weights and tail index k
    list(lwx, k)
  }

  if (!is.matrix(lw))
    lw <- as.matrix(lw)
  K <- ncol(lw)
  K2 <- 2 * K
  if (.Platform$OS.type != "windows") {
    vgis_out <- mclapply(X = 1:K, FUN = .vgis, mc.cores = cores)
  } else {
    cl <- makePSOCKcluster(cores)
    on.exit(stopCluster(cl))
    vgis_out <- parLapply(cl, X = 1:K, fun = .vgis)
  }
  # extract and return modified log weights and gPd shape param k estimates
  ux <- unlist(vgis_out, recursive = FALSE, use.names = FALSE)
  lw_smooth <- cbind_list(ux[seq(1, K2, 2)])
  pareto_k <- do.call(c, ux[seq(2, K2, 2)])
  nlist(lw_smooth, pareto_k)
}
