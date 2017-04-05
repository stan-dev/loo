#' Pareto smoothed importance sampling (PSIS)
#'
#' Implementation of Pareto smoothed importance sampling, a method for
#' stabilizing importance weights. For full details about the algorithm see
#' Vehtari, Gelman and Gabry (2016a, 2016b). For diagnostics see the
#' \link{pareto-k-diagnostic} page.
#'
#' @export
#' @param lw A matrix or vector of log weights. For computing LOO, \code{lw =
#'   -log_lik}, the \emph{negative} of an \eqn{S} (simulations) by \eqn{N} (data
#'   points) pointwise log-likelihood matrix.
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
#'   and \code{pareto_k} (estimated generalized
#'   Pareto \link[=pareto-k-diagnostic]{shape parameter(s) k}).
#'
#' @inheritSection loo-package PSIS-LOO
#'
#' @seealso \code{\link{pareto-k-diagnostic}} for PSIS diagnostics.
#'
#' @template loo-and-psis-references
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
    cutoff <- lw_cutpoint(x, wcp, MIN_CUTOFF)
    above_cut <- x > cutoff
    x_body <- x[!above_cut]
    x_tail <- x[above_cut]

    tail_len <- length(x_tail)
    if (tail_len < MIN_TAIL_LENGTH || all(x_tail == x_tail[1])) {
      if (all(x_tail == x_tail[1]))
        warning(
          "All tail values are the same. ",
          "Weights are truncated but not smoothed.",
          call. = FALSE
        )
      else if (tail_len < MIN_TAIL_LENGTH)
        warning(
          "Too few tail samples to fit generalized Pareto distribution.\n",
          "Weights are truncated but not smoothed.",
          call. = FALSE
        )

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
    # truncate (if wtrunc > 0) and renormalize,
    # return log weights and pareto k
    lw_new <- lw_normalize(lw_truncate(x_new, wtrunc))
    nlist(lw_new, k)
  }

  .psis_loop <- function(i) {
    if (LL_FUN) {
      ll_i <- llfun(i = i,
                    data = llargs$data[i,, drop=FALSE],
                    draws = llargs$draws)
      lw_i <- -1 * ll_i
    } else {
      lw_i <- lw[, i]
      ll_i <- -1 * lw_i
    }
    psis <- .psis(lw_i)
    if (FROM_LOO)
      nlist(lse = logSumExp(ll_i + psis$lw_new), k = psis$k)
    else
      psis
  }

  # minimal cutoff value. there must be at least 5 log-weights larger than this
  # in order to fit the gPd to the tail
  MIN_CUTOFF <- -700
  MIN_TAIL_LENGTH <- 5
  dots <- list(...)
  FROM_LOO <- if ("COMPUTE_LOOS" %in% names(dots))
    dots$COMPUTE_LOOS else FALSE

  if (!missing(lw)) {
    if (!is.matrix(lw))
      lw <- as.matrix(lw)
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
      # nocov start
      cl <- makePSOCKcluster(cores)
      on.exit(stopCluster(cl))
      out <- parLapply(cl, X = 1:N, fun = .psis_loop)
      # nocov end
    }
  }

  pareto_k <- vapply(out, "[[", 2L, FUN.VALUE = numeric(1))
  psislw_warnings(pareto_k)

  if (FROM_LOO) {
    loos <- vapply(out, "[[", 1L, FUN.VALUE = numeric(1))
    nlist(loos, pareto_k)
  } else {
    funval <- if (LL_FUN) llargs$S else nrow(lw)
    lw_smooth = vapply(out, "[[", 1L, FUN.VALUE = numeric(funval))
    nlist(lw_smooth, pareto_k)
  }
}


# internal ----------------------------------------------------------------
lw_cutpoint <- function(y, wcp, min_cut) {
  if (min_cut < log(.Machine$double.xmin))
    min_cut <- -700

  cp <- quantile(y, 1 - wcp, names = FALSE)
  max(cp, min_cut)
}

lw_truncate <- function(y, wtrunc) {
  if (wtrunc == 0)
    return(y)

  logS <- log(length(y))
  lwtrunc <- wtrunc * logS - logS + logSumExp(y)
  y[y > lwtrunc] <- lwtrunc
  y
}

#' @importFrom matrixStats logSumExp
lw_normalize <- function(y) {
  y - logSumExp(y)
}

# inverse-CDF of generalized Pareto distribution (formula from Wikipedia)
qgpd <- function(p, xi = 1, mu = 0, sigma = 1, lower.tail = TRUE) {
  if (is.nan(sigma) || sigma <= 0)
    return(rep(NaN, length(p)))
  if (!lower.tail)
    p <- 1 - p

  mu + sigma * ((1 - p)^(-xi) - 1) / xi
}


# warnings about pareto k values ------------------------------------------
psislw_warnings <- function(k) {
  if (any(k > 0.7)) {
    .warn(
      "Some Pareto k diagnostic values are too high. ",
      .k_help()
    )
  } else if (any(k > 0.5)) {
    .warn(
      "Some Pareto k diagnostic values are slightly high. ",
      .k_help()
    )
  }
}
