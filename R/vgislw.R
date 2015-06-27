#' Very good importance sampling (VGIS)
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

  # minimal cutoff value. there must be at least 5 log-weights larger than this
  # in order to fit the generalized Pareto distribution (gPd) to the tail
  MIN_CUTOFF <- -700
  MIN_TAIL_LENGTH <- 5

  .vgis <- function(i) {
    x <- lw[, i]
    S <- length(x)
    # split into body and right tail
    cutoff <- quantile(x, 1 - wcp, names = FALSE)
    cutoff <- max(cutoff, MIN_CUTOFF)
    x_cut <- x > cutoff
    xbody <- x[!x_cut]
    xtail <- x[x_cut]
    tail_len <- length(xtail)
    if (tail_len < MIN_TAIL_LENGTH) {
      # too few tail samples to fit gPd
      xnew <- x
      k <- Inf
    } else {
      # store order of tail samples
      tail_ord <- order(xtail)
      # fit gPd to the right tail samples
      xtail <- pmax(xtail, max(xtail) - thresh)
      exp_cutoff <- exp(cutoff)
      fit <- gpdfit(exp(xtail) - exp_cutoff)
      k <- min(fit$k, kmax)
      # compute order statistics for the fit
      qq <- qgpd(seq_min_half(tail_len)/tail_len, xi = k, beta = fit$sigma)
      qq <- qq + exp_cutoff
      # remap back to the original order
      smoothed_tail <- rep.int(0, tail_len)
      smoothed_tail[tail_ord] <- log(qq)
      # join body and gPd smoothed tail
      xnew <- x
      xnew[!x_cut] <- xbody
      xnew[x_cut] <- smoothed_tail
    }
    if (wtrunc > 0) {
      # truncate
      logS <- log(S)
      lwtrunc <- wtrunc * logS - logS + logSumExp(xnew)
      xnew[xnew > lwtrunc] <- lwtrunc
    }
    # renormalize weights
    lwx <- xnew - logSumExp(xnew)
    # return log weights and tail index k
    list(lwx, k)
  }

  stopifnot(is.numeric(lw))
  if (!is.matrix(lw))
    lw <- as.matrix(lw)
  K <- ncol(lw)
  if (.Platform$OS.type != "windows") {
    vgis_out <- mclapply(X = 1:K, FUN = .vgis, mc.cores = cores)
  } else {
    cl <- makePSOCKcluster(cores)
    on.exit(stopCluster(cl))
    vgis_out <- parLapply(cl, X = 1:K, fun = .vgis)
  }

  # Extract and return modified log weights and gPd shape param k estimates.
  # vgis_out is a list of length N=ncol(lw). each of the N elements of vgis_out
  # is itself a list of length 2. the first element is a vector of length
  # S=nrow(lw) containing the modified log weights and the second element is the
  # estimate of the pareto shape parameter k.
  ux <- unlist(vgis_out, recursive = FALSE, use.names = FALSE)
  # ux is now a list of length 2*N. the odd elements contain the modified log
  # weights and the even elements contain the pareto k estimates
  lw_smooth <- cbind_list(ux[nodds(K)])
  pareto_k <- unlist(ux[nevens(K)])
  nlist(lw_smooth, pareto_k)
}
