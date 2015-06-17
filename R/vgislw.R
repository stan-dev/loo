#' Very good importance sampling, log weights
#' @export
#' @param lw log weights.
#' @param wcp percentage of samples used for the generalized Pareto fit estimate.
#' @param wtrunc for truncating very large weights to n^\code{wtrunc}. No
#' trunction if \code{wtrunc} is \code{0}.
#' @param cores number of cores to use for parallelization.
#' @return a list with modified log weights and tail indices.
#'
vgislw <- function(lw, wcp = 20, wtrunc = 3/4, cores = parallel::detectCores()) {
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
      lwtrunc <- wtrunc * log(n) - log(n) + matrixStats::logSumExp(qx)
      qx[qx > lwtrunc] <- lwtrunc
    }
    # renormalize weights
    lwx <- qx - matrixStats::logSumExp(qx)
    # return log weights and tail index k
    list(lwx, fit$k)
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
