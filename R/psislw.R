#' Pareto smoothed importance sampling (PSIS)
#'
#' @export
#' @param lw A matrix or vector of log weights. For computing LOO, \code{lw =
#'   -log_lik} (see \code{\link{extract_log_lik}}) and is an \eqn{S} by \eqn{N}
#'   matrix where \eqn{S} is the number of simulations and \eqn{N} is the number
#'   of data points. (If \code{lw} is a vector it will be coerced to a
#'   one-column matrix.)
#' @param wcp The proportion of importance weights to use for the generalized
#'   Pareto fit. The \code{100*wcp}\% largest weights are used as the sample
#'   from which to estimate the parameters of the generalized Pareto
#'   distribution.
#' @param wtrunc For truncating very large weights to \eqn{S}^\code{wtrunc}. Set
#'   to zero for no truncation.
#' @param cores The number of cores to use for parallelization. This can be set
#'   for an entire R session by \code{options(loo.cores = NUMBER)}. The default
#'   is \code{\link[parallel]{detectCores}}().
#' @param llfun,llargs See \code{\link{loo.function}}.
#' @param ... Ignored when \code{psislw} is called directly. The \code{...} is
#'   only used internally when \code{psislw} is called by the \code{\link{loo}}
#'   function.
#'
#' @return A named list with components \code{lw_smooth} (modified log weights)
#'   and \code{pareto_k} (estimated generalized Pareto shape parameter(s)
#'   \eqn{k}).
#'
#' @details See the 'PSIS-LOO' section in \code{\link{loo-package}}.
#'
#' @template internal-function-note
#' @template loo-paper-reference
#'
#' @importFrom matrixStats logSumExp
#' @importFrom parallel mclapply makePSOCKcluster stopCluster parLapply
#'
psislw <- function(lw, wcp = 0.2, wtrunc = 3/4,
                   cores = getOption("loo.cores", parallel::detectCores()),
                   llfun = NULL, llargs = NULL,
                   ...) {
  .psis <- function(lw_i) {
    x <- lw_i - max(lw_i)
    # split into body and right tail
    cutoff <- lw_cutpoint(x, wcp, MIN_CUTOFF)
    above_cut <- x > cutoff
    x_body <- x[!above_cut]
    x_tail <- x[above_cut]
    tail_len <- length(x_tail)
    if (tail_len < MIN_TAIL_LENGTH) {
      # too few tail samples to fit gPd
      warning("Too few tail samples to fit generalized Pareto distribution.\n",
              "Weights are truncated and normalized but not smoothed.",
              call. = FALSE)
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
    lw_new <- lw_normalize(lw_truncate(x_new, wtrunc))
    nlist(lw_new, k)
  }

  .psis_loop <- function(i) {
    if (LL_FUN) {
      ll_i <- llfun(i = i, data = llargs$data[i,, drop=FALSE],
                    draws = llargs$draws)
      lw_i <- -1 * ll_i
    } else {
      lw_i <- lw[, i]
      ll_i <- -1 * lw_i
    }
    psis <- .psis(lw_i)
    if (!FROM_LOO) psis
    else nlist(lse = logSumExp(ll_i + psis$lw_new), k = psis$k)
  }

  # minimal cutoff value. there must be at least 5 log-weights larger than this
  # in order to fit the gPd to the tail
  MIN_CUTOFF <- -700
  MIN_TAIL_LENGTH <- 5
  dots <- list(...)
  FROM_LOO <- if ("COMPUTE_LOOS" %in% names(dots))
    dots$COMPUTE_LOOS else FALSE

  if (!missing(lw)) {
    if (!is.matrix(lw)) lw <- as.matrix(lw)
    N <- ncol(lw)
    LL_FUN <- FALSE
  } else {
    if (is.null(llfun) || is.null(llargs))
      stop("Either 'lw' or 'llfun' and 'llargs' must be specified.")
    N <- llargs$N
    LL_FUN <- TRUE
  }
  if (cores == 1) {
    # don't call functions from parallel package if cores=1
    out <- lapply(X = 1:N, FUN = .psis_loop)
  } else {
    # parallelize
    if (.Platform$OS.type != "windows") {
      out <- mclapply(X = 1:N, FUN = .psis_loop, mc.cores = cores)
    } else {
      cl <- makePSOCKcluster(cores)
      on.exit(stopCluster(cl))
      out <- parLapply(cl, X = 1:N, fun = .psis_loop)
    }
  }
  pareto_k <- vapply(out, "[[", 2L, FUN.VALUE = numeric(1))
  if (FROM_LOO) {
    nlist(loos = vapply(out, "[[", 1L, FUN.VALUE = numeric(1)),
          pareto_k)
  } else {
    funval <- if (LL_FUN) llargs$S else nrow(lw)
    nlist(lw_smooth = vapply(out, "[[", 1L, FUN.VALUE = numeric(funval)),
          pareto_k)
  }
}
