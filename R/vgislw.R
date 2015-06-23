#' @export
#' @rdname vgisloo
#' @param lw an \eqn{S} by \eqn{N} matrix of log weights (\code{-log_lik} for
#'   computing LOO).
#' @param wcp the percentage of samples used for the generalized Pareto fit
#'   estimate.
#' @param wtrunc for truncating very large weights to \eqn{N}^\code{wtrunc}. Set
#'   to zero for no truncation.
#' @param fix_value if the largest value in any column of \code{lw} is greater
#'   than the second largest by at least \code{fix_value}, then the second
#'   largest value will be set equal to the largest. Set to zero to skip this
#'   step.
#' @param cores the number of cores to use for parallelization.
#'
vgislw <- function(lw, wcp = 20, wtrunc = 3/4, fix_value = 100,
                   cores = parallel::detectCores()) {
  .loop_fn <- function(i) {
    x <- lw[, i]
    # divide log weights into body and right tail
    n <- length(x)
    cutoff <- quantile(x, 1 - wcp / 100, names = FALSE)
    x_cut <- x > cutoff
    x1 <- x[!x_cut]
    x2 <- x[x_cut]
    n2 <- length(x2)
    # store order of tail samples
    x2_order <- order(x2)
    # fit generalized Pareto distribution to the right tail samples
    exp_cutoff <- exp(cutoff)
    fit <- gpdfit(exp(x2) - exp_cutoff)
    # compute ordered statistic for the fit
    qq <- qgpd(seq_min_half(n2) / n2, xi = fit$k, beta = fit$sigma) + exp_cutoff
    # remap back to the original order
    slq <- rep.int(0, n2)
    slq[x2_order] <- log(qq)
    # join body and GPD smoothed tail
    qx <- x
    qx[!x_cut] <- x1
    qx[x_cut] <- slq
    if (wtrunc > 0) {
      # truncate too large weights
      logn <- log(n)
      lwtrunc <- wtrunc * logn - logn + matrixStats::logSumExp(qx)
      qx[qx > lwtrunc] <- lwtrunc
    }
    # renormalize weights
    lwx <- qx - matrixStats::logSumExp(qx)
    # return log weights and tail index k
    list(lwx, fit$k)
  }

  if (fix_value > 0) {
    lw <- fix_large_diffs(lw, fix_value = fix_value)
  }

  K <- ncol(lw)
  K2 <- 2 * K
  if (.Platform$OS.type != "windows") {
    out <- parallel::mclapply(1:K, .loop_fn, mc.cores = cores)
  } else {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl))
    out <- parallel::parLapply(cl, X = 1:K, fun = .loop_fn)
  }
  ux <- unlist(out, recursive = FALSE, use.names = FALSE)
  lw <- do.call(cbind, ux[seq(1, K2, 2)])
  kss <- do.call(c, ux[seq(2, K2, 2)])

  # return log weights and tail indices k
  nlist(lw, k = kss)
}
