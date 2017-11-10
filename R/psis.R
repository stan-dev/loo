#' Pareto smoothed importance sampling (PSIS)
#'
#' Implementation of Pareto smoothed importance sampling, a method for
#' stabilizing importance weights. For full details about the algorithm see
#' Vehtari, Gelman and Gabry (2017a, 2017b). For diagnostics see the
#' \link{pareto-k-diagnostic} page.
#'
#' @export
#' @param x An array, matrix, or vector of log-weights (for LOO these are
#'   negative log-likelihood values). See the \strong{Methods (by class)}
#'   section below for a detailed description of how to specify the inputs for
#'   each method.
#' @param ... Arguments passed on to the various methods.
#' @param wtrunc Should weights be truncated? Set to \code{FALSE} or \code{0}
#'   no truncation. As of version \code{2.0.0}, the default is to truncate
#'   the smoothed weights at the (column) maximum(s) of the raw weights.
#' @param cores The number of cores to use for parallelization. The default for
#'   an entire R session can be set with \code{options(loo.cores = NUMBER)}. As
#'   of version \code{2.0.0} the \strong{default is now 1 core}, but we
#'   recommend using as many (or close to as many) cores as possible.
#'
#' @return The \code{psis} methods return an object of class \code{"psis"},
#'   which is a named list with the following components:
#'
#' \describe{
#'   \item{\code{log_weights}}{
#'     Vector or matrix of smoothed but \emph{unnormalized} log weights. To get
#'     normalized weights use the \code{weights} method provided for objects
#'     of class \code{"psis"}.
#'   }
#'   \item{\code{pareto_k}}{
#'     Vector of estimated \link[=pareto-k-diagnostic]{shape parameter(s) k} of
#'     the generalized Pareto distribution.
#'   }
#'   \item{\code{n_eff}}{
#'     Vector of estimated PSIS effective sample size(s).
#'   }
#' }
#'
#' Objects of class \code{"psis"} also have the following
#' \code{\link{attributes}}:
#' \describe{
#'   \item{\code{norm_const_log}}{
#'     Vector of precomputed values of \code{colLogSumExps(log_weights)} that are
#'     used internally by the \code{weights} method to normalize the log weights.
#'   }
#'   \item{\code{tail_len}}{
#'     Vector of tail lengths used for fitting the generalized Pareto
#'     distribution.
#'   }
#'   \item{\code{rel_n_eff}}{
#'     Vector of relative effective sample size estimates of the exponentiated
#'     log-likelihood.
#'   }
#'   \item{\code{dims}}{
#'     Integer vector of length 2 containing \code{S} (posterior sample size)
#'     and \code{N} (number of observations).
#'   }
#' }
#'
#' @seealso \code{\link{pareto-k-diagnostic}} for PSIS diagnostics.
#'
#' @template loo-and-psis-references
#'
psis <- function(x, ...) UseMethod("psis")

#' @export
#' @templateVar fn psis
#' @template array
#'
psis.array <-
  function(x,
           ...,
           cores = getOption("loo.cores", 1),
           wtrunc = TRUE) {
    stopifnot(length(dim(x)) == 3)
    lw <- validate_ll(x)
    rel_n_eff <- relative_n_eff(exp(-lw))

    lw <- llarray_to_matrix(lw)
    do_psis(
      lw = lw,
      rel_n_eff = rel_n_eff,
      wtrunc = wtrunc,
      cores = cores
    )
  }

#' @export
#' @templateVar fn psis
#' @template matrix
#' @template chain_id
#'
psis.matrix <-
  function(x,
           chain_id,
           ...,
           cores = getOption("loo.cores", 1),
           wtrunc = TRUE) {
    lw <- validate_ll(x)
    rel_n_eff <- relative_n_eff(exp(-lw), chain_id)
    do_psis(
      lw = lw,
      rel_n_eff = rel_n_eff,
      wtrunc = wtrunc,
      cores = cores
    )
  }

#' @export
#' @describeIn psis
#' A vector of length \eqn{S} (posterior sample size).
#'
psis.default <-
  function(x,
           chain_id,
           ...,
           wtrunc = TRUE) {
    stopifnot(is.null(dim(x)) || length(dim(x)) == 1)
    dim(x) <- c(length(x), 1)
    psis.matrix(x = x,
                chain_id = chain_id,
                wtrunc = wtrunc,
                cores = 1)
  }

#' @rdname psis
#' @export
#' @method weights psis
#' @param object For the \code{weights} method, an object returned by
#'   \code{psis} (a list with class \code{"psis"}).
#' @param log For the \code{weights} method, should the weights be returned on
#'   the log scale? Defaults to \code{TRUE}.
#' @param normalize For the \code{weights} method, should the weights be
#'   normalized? Defaults to \code{TRUE}.
#'
#' @return The \code{weights} method returns an object with the same dimensions
#'   as the \code{log_weights} component of the \code{"psis"} object. The
#'   \code{normalize} and \code{log} arguments control whether the returned
#'   weights are normalized and whether or not to return them on the log scale.
#'
weights.psis <-
  function(object,
           ...,
           log = TRUE,
           normalize = TRUE) {
    out <- object[["log_weights"]]  # smoothed but unnormalized log weights
    const <- attr(object, "norm_const_log")  # precomputed colLogSumExp(log_weights)
    if (normalize)
      out <- sweep(out, 2, const)
    if (!log)
      out <- exp(out)

    return(out)
  }


# internal ----------------------------------------------------------------

# Do PSIS given matrix of log weights
#
# @param lw Matrix of log weights (-loglik)
# @param rel_n_eff Vector of relative effective sample sizes
# @param wtrunc User's scalar wtrunc argument
# @param cores User's integer cores argument
# @return A list with class "psis" and structure described in the main doc at
#   the top of this file.
#
do_psis <- function(lw, rel_n_eff = 1, wtrunc = TRUE, cores) {
  stopifnot(length(rel_n_eff) == ncol(lw),
            cores == as.integer(cores),
            wtrunc >= 0)
  N <- ncol(lw)
  S <- nrow(lw)
  tail_len <- n_pareto(rel_n_eff, S)

  if (cores == 1) {
    lw_list <- lapply(seq_len(N), function(i)
      do_psis_i(lw[, i], tail_len[i]))
  } else {
    if (.Platform$OS.type != "windows") {
      lw_list <- parallel::mclapply(
        X = seq_len(N),
        FUN = function(i)
          do_psis_i(lw[, i], tail_len[i]),
        mc.cores = cores
      )
    } else {
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(parallel::stopCluster(cl))
      lw_list <-
        parallel::parLapply(
          cl = cl,
          X = seq_len(N),
          fun = function(i)
            do_psis_i(lw[, i], tail_len[i])
        )
    }
  }

  pareto_k <- psis_apply(lw_list, "pareto_k")
  throw_psis_warnings(pareto_k)

  log_weights <- psis_apply(lw_list, "log_weights", fun_val = numeric(S))
  if (wtrunc) {
    # truncate at the max of raw wts
    max_wts <- apply(lw, 2, max)
    log_weights <- sapply(1:N, function(n) {
      truncate_lw(log_weights[, n], upper = max_wts[n])
    })
  }

  psis_object(
    unnormalized_log_weights = log_weights,
    pareto_k = pareto_k,
    tail_len = tail_len,
    rel_n_eff = rel_n_eff,
    wtrunc = wtrunc
  )
}

#' Extract named components from each list in a list of lists
#'
#' @noRd
#' @param x List of lists.
#' @param item String naming the component or attribute to pull out of each list (or list-like object).
#' @param fun,fun.val passed to vapply's FUN and FUN.VALUE.
#' @return Numeric vector or matrix.
#'
psis_apply <- function(x, item, fun = c("[[", "attr"), fun_val = numeric(1)) {
  stopifnot(is.list(x))
  vapply(x, FUN = match.arg(fun), FUN.VALUE = fun_val, item)
}

#' Structure the object returned by the psis methods
#'
#' @noRd
#' @param unnormalized_log_weights Smoothed and possibly truncated log weights,
#'   but unnormalized.
#' @param pareto_k Vector of GPD k estimates.
#' @param tail_len Vector of tail lengths used to fit GPD.
#' @param rel_n_eff Vector of relative MCMC n_eff for exp(log lik)
#' @param wtrunc User's wtrunc argument.
#' @return A list of class "psis" with structure described in the main doc at the
#'   top of this file.
#'
psis_object <-
  function(unnormalized_log_weights,
           pareto_k,
           tail_len,
           rel_n_eff,
           wtrunc) {
    stopifnot(is.matrix(unnormalized_log_weights))

    lwdim <- setNames(dim(unnormalized_log_weights), c("S", "N"))
    norm_const_log <- colLogSumExps(unnormalized_log_weights)

    out <- structure(
      list(
        log_weights = unnormalized_log_weights,
        pareto_k = pareto_k
      ),
      # attributes
      norm_const_log = norm_const_log,
      tail_len = tail_len,
      rel_n_eff = rel_n_eff,
      wtrunc = wtrunc,
      dims = lwdim,
      class = c("psis", "list")
    )

    # need weights (normalized and not log weights) to compute psis_n_eff
    w <- weights(out, normalize = TRUE, log = FALSE)
    out[["n_eff"]] <- psis_n_eff(w, rel_n_eff)
    return(out)
  }


#' PSIS (without truncation and normalization) on a single vector
#'
#' @noRd
#' @param lw_i A vector of log weights (negative log likelihood).
#' @param tail_len_i An integer tail length.
#' @return list containing lw (vector of untruncated and unnormalized log
#'   weights), and pareto_k ('khat' estimate).
#'
do_psis_i <- function(lw_i, tail_len_i) {
  S <- length(lw_i)
  lw_i <- lw_i - max(lw_i)
  khat <- Inf

  if (!enough_tail_samples(tail_len_i)) {
    warning(
      "Too few tail samples to fit generalized Pareto distribution. ",
      "Weights will be truncated (if wtrunc > 0) but not smoothed.",
      call. = FALSE
    )
  } else {
    tail_ids <- seq(S - tail_len_i + 1, S)
    ord <- sort.int(lw_i, index.return = TRUE)
    lw_tail <- ord$x[tail_ids]
    if (all(lw_tail == lw_tail[1])) {
      warning(
        "All tail values are the same. ",
        "Weights will be truncated (if wtrunc > 0) but not smoothed.",
        call. = FALSE
      )
    } else {
      smoothed <- psis_smooth_tail(lw_tail)
      ord$x[tail_ids] <- smoothed$tail
      khat <- smoothed$k
      lw_i <- ord$x[ord$ix]
    }
  }
  list(log_weights = lw_i, pareto_k = khat)
}


#' PSIS tail smoothing for a single vector
#'
#' @noRd
#' @param x Vector of tail elements already sorted in ascending order.
#' @return List with components 'tail' (vector same size as x) and 'k' (scalar
#'   shape parameter estimate).
#' The 'tail' component contains the logs of the order statistics of the GPD.
#'
psis_smooth_tail <- function(x) {
  len <- length(x)
  exp_cutoff <- exp(x[1])
  fit <- gpdfit(exp(x) - exp_cutoff)
  k <- fit$k
  sigma <- fit$sigma

  p <- (seq_len(len) - 0.5) / len
  qq <- qgpd(p, k, sigma) + exp_cutoff
  list(tail = log(qq), k = k)
}


# Calculate tail lengths to use for fitting the GPD
#
# The number of weights (i.e., tail length) used to fit the generalized Pareto
# distribution is now decreasing with the number of posterior draws S, and is
# also adjusted based on the relative MCMC neff for exp(log_lik). This will
# answer the questions about the asymptotic properties, works better for thick
# tailed proposal distributions, and is adjusted based on dependent Markov chain
# samples. Specifically, the tail length is now 3*sqrt(S)/rel_neff but capped at
# 20% of the total number of weights.
#
# @param rel_n_eff A N-vector of relative MCMC effective sample sizes of
#   exp(log-lik matrix).
# @param S The (integer) size of posterior sample.
# @return An N-vector of tail lengths.
#
n_pareto <- function(rel_n_eff, S) {
  n <- pmin(0.2 * S, 3 * sqrt(S / rel_n_eff))
  ceiling(n)
}

# Check for enough tail samples to fit GPD
#
# @param tail_len Integer tail length.
# @param min_len The minimum allowed tail length.
# @return TRUE or FALSE
#
enough_tail_samples <- function(tail_len, min_len = 5) {
  tail_len >= min_len
}


# Truncate a vector of (log) weights
#
# @param x A vector of log weights.
# @param logS Precomputed scalar log(length(x)).
# @param wtrunc User's scalar wtrunc parameter.
# @return x, but modified so all values above the truncation value are set equal
#   to the truncation value.
#
truncate_lw <- function(x, upper) {
  x[x > upper] <- upper
  return(x)
}
# truncate_lw_old <- function(x, logS, wtrunc) {
#   val <- wtrunc * logS - logS + logSumExp(x)
#   x[x > val] <- val
#   return(x)
# }


#' Throw warnings about pareto k estimates
#'
#' @noRd
#' @param k A vector of pareto k estimates.
#' @param high The value at which to warn about slighly high estimates.
#' @param too_high The value at which to warn about very high estimates.
#' @return Nothing, just possibly throws warnings.
#'
throw_psis_warnings <- function(k, high = 0.5, too_high = 0.7) {
  if (any(k > too_high)) {
    .warn("Some Pareto k diagnostic values are too high. ",
          .k_help())
  } else if (any(k > high)) {
    .warn("Some Pareto k diagnostic values are slightly high. ",
          .k_help())
  }
}

