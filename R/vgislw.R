#' Very good importance sampling, log weights
#' @export
#' @param lw log weights.
#' @param wcp percentage of samples used for the generalized Pareto fit estimate.
#' @param wtrunc for truncating very large weights to \code{n^wtrunc}. No
#' trunction if \code{wtrunc=0}.
#' @param cores number of cores to use for parallelization.
#' @return a list with modified log waits and tail indices.
#'
vgislw <- function(lw, wcp = 20, wtrunc = 3/4, cores = parallel::detectCores()) {
  .loop_fn <- function(i) {
    x <- lw[, i]
    # divide log weights into body and right tail
    n <- length(x)
    cutoff <- quantile(x, 1 - wcp / 100, names = FALSE)
    x_gt_cut <- x > cutoff
    x1 <- x[!x_gt_cut]
    x2 <- x[x_gt_cut]
    n2 <- length(x2)
    # store order of tail samples
    x2si <- order(x2)
    # fit generalized Pareto distribution to the right tail samples
    fit <- gpdfit(exp(x2) - exp(cutoff))
    k <- fit$k
    sigma <- fit$sigma
    # compute ordered statistic for the fit
    qq <- qgpd(seq_min_half(n2) / n2, xi = k, beta = sigma) + exp(cutoff)
    # remap back to the original order
    slq <- rep.int(0, n2)
    slq[x2si] <- log(qq)
    # join body and GPD smoothed tail
    qx <- x
    qx[!x_gt_cut] <- x1
    qx[x_gt_cut] <- slq
    if (wtrunc > 0) {
      # truncate too large weights
      lwtrunc <- wtrunc * log(n) - log(n) + sumlogs(qx)
      qx[qx > lwtrunc] <- lwtrunc
    }

    # renormalize weights
    lwx <- qx - sumlogs(qx)

    # return log weights and tail index k
    list(lwx, k)
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
