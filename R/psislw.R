#' Pareto smoothed importance sampling (PSIS)
#'
#' @export
#' @param lw a matrix or vector of log weights. For computing LOO \code{lw =
#'   -log_lik} (see \code{\link{extract_log_lik}}) and is an \eqn{S} by \eqn{N}
#'   matrix where \eqn{S} is the number of simulations and \eqn{N} is the number
#'   of data points. (If \code{lw} is a vector it will be coerced to a
#'   one-column matrix.)
#' @param wcp the proportion of importance weights to use for the generalized
#'   Pareto fit. The \code{100*wcp}\% largest weights are used as the sample
#'   from which to estimate the parameters of the generalized Pareto
#'   distribution.
#' @param wtrunc for truncating very large weights to \eqn{S}^\code{wtrunc}. Set
#'   to zero for no truncation.
#' @param cores the number of cores to use for parallelization.
#'
#' @return A named with list with components \code{lw_smooth} (modified log
#'   weights) and \code{pareto_k} (estimated generalized Pareto shape parameters
#'   \eqn{k}).
#'
#' @details See the 'PSIS-LOO' section in \code{\link{loo-package}}.
#'
#' @note This function is primarily intended for internal use, but is exported
#'   so that users can call it directly for other purposes. Users simply
#'   wishing to compute LOO should use the \code{\link{loo}} function.
#'
#' @seealso \code{\link{loo-package}}, \code{\link{loo}}
#'
#' @importFrom matrixStats logSumExp
#' @importFrom parallel mclapply makePSOCKcluster stopCluster parLapply
#'
psislw <- function(lw, wcp = 0.2, wtrunc = 3/4,
                   cores = parallel::detectCores()) {
  .psis <- function(n) {
    x <- lw[, n]
    x <- x - max(x)
    # split into body and right tail
    cutoff <- lw_cutpoint(x, wcp, MIN_CUTOFF)
    above_cut <- x > cutoff
    x_body <- x[!above_cut]
    x_tail <- x[above_cut]
    tail_len <- length(x_tail)
    if (tail_len < MIN_TAIL_LENGTH) {
      # too few tail samples to fit gPd
      x_new <- x
      k <- Inf
    } else {
      # store order of tail samples, fit gPd to the right tail samples, compute
      # order statistics for the fit, remap back to the original order, join
      # body and gPd smoothed tail
      tail_ord <- order(x_tail)
      exp_cutoff <- exp(cutoff)
      fit <- gpdfit(exp(x_tail) - exp_cutoff)
      k <- fit$k
      prb <- (seq_len(tail_len) - 0.5) / tail_len
      qq <- qgpd(p = prb, xi = k, sigma = fit$sigma) + exp_cutoff
      smoothed_tail <- rep.int(0, tail_len)
      smoothed_tail[tail_ord] <- log(qq)
      x_new <- x
      x_new[!above_cut] <- x_body
      x_new[above_cut] <- smoothed_tail
    }
    # truncate (if wtrunc > 0) and renormalize, return log weights and pareto k
    lw_new <- lw_truncate(x_new, wtrunc)
    lw_new <- lw_normalize(lw_new)
    nlist(lw_new, k)
  }

  # minimal cutoff value. there must be at least 5 log-weights larger than this
  # in order to fit the gPd to the tail
  MIN_CUTOFF <- -700
  MIN_TAIL_LENGTH <- 5

  if (!is.matrix(lw))
    lw <- as.matrix(lw)
  N <- ncol(lw)
  if (.Platform$OS.type != "windows") {
    psis <- mclapply(X = 1:N, FUN = .psis, mc.cores = cores)
  } else {
    cl <- makePSOCKcluster(cores)
    on.exit(stopCluster(cl))
    psis <- parLapply(cl, X = 1:N, fun = .psis)
  }
  .psis_out(psis)
}
