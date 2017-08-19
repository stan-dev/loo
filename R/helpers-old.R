
logColMeansExp_llfun <- function(fun, args) {
  # should be more stable than log(colMeans(exp(x)))
  logS <- log(args$S)
  colLSEs <- vapply(seq_len(args$N), FUN = function(i) {
    logSumExp(fun(i = i, data = args$data[i,,drop=FALSE], draws = args$draws))
  }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
  colLSEs - logS
}

old_pointwise_loo <- function(psis, log_lik, llfun = NULL, llargs = NULL) {
  if (!missing(log_lik)) {
    lpd <- logColMeansExp(log_lik)
  } else {
    if (is.null(llfun) || is.null(llargs))
      stop("Either 'log_lik' or 'llfun' and 'llargs' must be specified.",
           call. = FALSE)
    lpd <- logColMeansExp_llfun(llfun, llargs)
  }
  elpd_loo <- psis$loos
  p_loo <- lpd - elpd_loo
  looic <- -2 * elpd_loo
  pointwise <- nlist(elpd_loo, p_loo, looic)
  out <- totals(pointwise)
  nms <- names(pointwise)
  names(out) <- c(nms, paste0("se_", nms))
  out$pointwise <- do.call(cbind, pointwise)
  out$pareto_k <- psis$pareto_k
  out
}

# Compute estimates from pointwise vectors
#
# @param pointwise A list of vectors.
# @return A named list of estimates and standard errors.
#
totals <- function(pointwise) {
  N <- length(pointwise[[1L]])
  ests <- lapply(pointwise, sum)
  ses <- lapply(pointwise, function(x) sqrt(N * var(x)))
  names(ses) <- paste0("se_", names(ses))
  c(ests, ses)
}

