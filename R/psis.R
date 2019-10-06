#' Pareto smoothed importance sampling (PSIS)
#'
#' Implementation of Pareto smoothed importance sampling (PSIS), a method for
#' stabilizing importance ratios. The version of PSIS implemented here
#' corresponds to the algorithm presented in Vehtari, Gelman and Gabry (2017b).
#' For PSIS diagnostics see the [pareto-k-diagnostic] page.
#'
#' @export
#' @param log_ratios An array, matrix, or vector of importance ratios on the log
#'   scale (for PSIS-LOO these are *negative* log-likelihood values). See the
#'   **Methods (by class)** section below for a detailed description of how
#'   to specify the inputs for each method.
#' @param ... Arguments passed on to the various methods.
#' @template cores
#' @param r_eff Vector of relative effective sample size estimates containing
#'   one element per observation. The values provided should be the relative
#'   effective sample sizes of `1/exp(log_ratios)` (i.e., `1/ratios`).
#'   This is related to the relative efficiency of estimating the normalizing
#'   term in self-normalizing importance sampling. If `r_eff` is not
#'   provided then the reported PSIS effective sample sizes and Monte Carlo
#'   error estimates will be over-optimistic. See the [relative_eff()]
#'   helper function for computing `r_eff`. If using `psis` with
#'   draws of the `log_ratios` not obtained from MCMC then the warning
#'   message thrown when not specifying `r_eff` can be disabled by
#'   setting `r_eff` to `NA`.
#'
#' @return The `psis()` methods return an object of class `"psis"`,
#'   which is a named list with the following components:
#'
#' \describe{
#'   \item{`log_weights`}{
#'     Vector or matrix of smoothed (and truncated) but *unnormalized* log
#'     weights. To get normalized weights use the `weights` method provided
#'     for objects of class `"psis"`.
#'   }
#'  \item{`diagnostics`}{
#'    A named list containing two vectors:
#'    * `pareto_k`: Estimates of the shape parameter \eqn{k} of the
#'      generalized Pareto distribution. See the [pareto-k-diagnostic]
#'      page for details.
#'    * `n_eff`: PSIS effective sample size estimates.
#'  }
#' }
#'
#' Objects of class `"psis"` also have the following [attributes][attributes()]:
#' \describe{
#'   \item{`norm_const_log`}{
#'     Vector of precomputed values of `colLogSumExps(log_weights)` that are
#'     used internally by the `weights` method to normalize the log weights.
#'   }
#'   \item{`tail_len`}{
#'     Vector of tail lengths used for fitting the generalized Pareto distribution.
#'   }
#'   \item{`r_eff`}{
#'     If specified, the user's `r_eff` argument.
#'   }
#'   \item{`dims`}{
#'     Integer vector of length 2 containing `S` (posterior sample size)
#'     and `N` (number of observations).
#'   }
#' }
#'
#' @seealso
#' \itemize{
#' \item [loo()] for approximate LOO-CV using PSIS.
#' \item [pareto-k-diagnostic] for PSIS diagnostics.
#' }
#'
#' @template loo-and-psis-references
#'
#' @examples
#' log_ratios <- -1 * example_loglik_array()
#' r_eff <- relative_eff(exp(-log_ratios))
#' psis_result <- psis(log_ratios, r_eff = r_eff)
#' str(psis_result)
#' plot(psis_result)
#'
#' # extract smoothed weights
#' lw <- weights(psis_result) # default args are log=TRUE, normalize=TRUE
#' ulw <- weights(psis_result, normalize=FALSE) # unnormalized log-weights
#'
#' w <- weights(psis_result, log=FALSE) # normalized weights (not log-weights)
#' uw <- weights(psis_result, log=FALSE, normalize = FALSE) # unnormalized weights
#'
#'
psis <- function(log_ratios, ...) UseMethod("psis")


#' @export
#' @templateVar fn psis
#' @template array
#'
psis.array <-
  function(log_ratios, ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
    cores <- loo_cores(cores)
    stopifnot(length(dim(log_ratios)) == 3)
    log_ratios <- validate_ll(log_ratios)
    log_ratios <- llarray_to_matrix(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_psis(log_ratios, r_eff = r_eff, cores = cores)
  }


#' @export
#' @templateVar fn psis
#' @template matrix
#'
psis.matrix <-
  function(log_ratios,
           ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
    cores <- loo_cores(cores)
    log_ratios <- validate_ll(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_psis(log_ratios, r_eff = r_eff, cores = cores)
  }


#' @export
#' @templateVar fn psis
#' @template vector
#'
psis.default <-
  function(log_ratios, ..., r_eff = NULL) {
    stopifnot(is.null(dim(log_ratios)) || length(dim(log_ratios)) == 1)
    dim(log_ratios) <- c(length(log_ratios), 1)
    r_eff <- prepare_psis_r_eff(r_eff, len = 1)
    psis.matrix(log_ratios, r_eff = r_eff, cores = 1)
  }


#' @rdname psis
#' @export
#' @export weights.psis
#' @method weights psis
#' @param object,x An object returned by `psis()`.
#' @param log For the `weights()` method, should the weights be returned on
#'   the log scale? Defaults to `TRUE`.
#' @param normalize For the `weights()` method, should the weights be
#'   normalized? Defaults to `TRUE`.
#'
#' @return The `weights()` method returns an object with the same dimensions as
#'   the `log_weights` component of the `"psis"` object. The `normalize` and
#'   `log` arguments control whether the returned weights are normalized and
#'   whether or not to return them on the log scale.
#'
weights.psis <-
  function(object,
           ...,
           log = TRUE,
           normalize = TRUE) {
    out <- object[["log_weights"]] # smoothed but unnormalized log weights
    if (normalize) {
      out <- sweep(out,
                   MARGIN = 2,
                   STATS = attr(object, "norm_const_log"), # colLogSumExp(log_weights)
                   check.margin = FALSE)
    }
    if (!log) {
      out <- exp(out)
    }

    return(out)
  }


# Subset a psis object without breaking it
#
#' @rdname psis
#' @export
#' @export subset.psis
#' @method subset psis
#' @param subset For the `subset()` method, a vector indicating which
#'   observations (columns of weights) to keep. Can be a logical vector of
#'   length `ncol(x)` (for a psis object `x`) or a shorter integer vector
#'   containing only the indexes to keep.
#'
#' @return The `subset()` returns a `"psis"` object. It is the same as the input
#'   but without the contents corresponding to the unselected indexes.
#'
subset.psis <- function(x, subset, ...) {
  if (anyNA(subset)) {
    stop("NAs not allowed in subset.", call. = FALSE)
  }
  if (is.logical(subset) || all(subset %in% c(0,1))) {
    stopifnot(length(subset) == dim(x)[2])
    subset <- which(as.logical(subset))
  } else {
    stopifnot(length(subset) <= dim(x)[2],
              all(subset == as.integer(subset)))
    subset <- as.integer(subset)
  }

  x$log_weights <- x$log_weights[, subset, drop=FALSE]
  x$diagnostics$pareto_k <- x$diagnostics$pareto_k[subset]
  x$diagnostics$n_eff <- x$diagnostics$n_eff[subset]

  structure(
    .Data = x,
    class = class(x),
    dims = c(dim(x)[1], length(subset)),
    norm_const_log = attr(x, "norm_const_log")[subset],
    tail_len = attr(x, "tail_len")[subset],
    r_eff = attr(x, "r_eff")[subset],
    subset = subset
  )
}


#' @rdname psis
#' @export
dim.psis <- function(x) {
  attr(x, "dims")
}


#' @rdname psis
#' @export
is.psis <- function(x) {
  inherits(x, "psis") && is.list(x)
}



# internal ----------------------------------------------------------------

#' Structure the object returned by the psis methods
#'
#' @noRd
#' @param unnormalized_log_weights Smoothed and possibly truncated log weights,
#'   but unnormalized.
#' @param pareto_k Vector of GPD k estimates.
#' @param tail_len Vector of tail lengths used to fit GPD.
#' @param r_eff Vector of relative MCMC n_eff for `exp(log lik)`
#' @return A list of class `"psis"` with structure described in the main doc at
#'   the top of this file.
#'
psis_object <-
  function(unnormalized_log_weights,
           pareto_k,
           tail_len,
           r_eff) {
    stopifnot(is.matrix(unnormalized_log_weights))
    norm_const_log <- matrixStats::colLogSumExps(unnormalized_log_weights)
    out <- structure(
      list(
        log_weights = unnormalized_log_weights,
        diagnostics = list(pareto_k = pareto_k, n_eff = NULL)
      ),
      # attributes
      norm_const_log = norm_const_log,
      tail_len = tail_len,
      r_eff = r_eff,
      dims = dim(unnormalized_log_weights),
      class = c("psis", "list")
    )

    # need normalized weights (not on log scale) for psis_n_eff
    w <- weights.psis(out, normalize = TRUE, log = FALSE)
    out$diagnostics[["n_eff"]] <- psis_n_eff(w, r_eff)
    return(out)
  }


#' Do PSIS given matrix of log weights
#'
#' @noRd
#' @param lr Matrix of log ratios (`-loglik`)
#' @param r_eff Vector of relative effective sample sizes
#' @param cores User's integer `cores` argument
#' @return A list with class `"psis"` and structure described in the main doc at
#'   the top of this file.
#'
do_psis <- function(log_ratios, r_eff, cores) {
  stopifnot(cores == as.integer(cores))
  N <- ncol(log_ratios)
  S <- nrow(log_ratios)
  tail_len <- n_pareto(r_eff, S)
  throw_tail_length_warnings(tail_len)

  if (cores == 1) {
    lw_list <- lapply(seq_len(N), function(i)
      do_psis_i(log_ratios[, i], tail_len[i]))
  } else {
    if (.Platform$OS.type != "windows") {
      lw_list <- parallel::mclapply(
        X = seq_len(N),
        mc.cores = cores,
        FUN = function(i)
          do_psis_i(log_ratios[, i], tail_len[i])
      )
    } else {
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(parallel::stopCluster(cl))
      lw_list <-
        parallel::parLapply(
          cl = cl,
          X = seq_len(N),
          fun = function(i)
            do_psis_i(log_ratios[, i], tail_len[i])
        )
    }
  }

  log_weights <- psis_apply(lw_list, "log_weights", fun_val = numeric(S))
  pareto_k <- psis_apply(lw_list, "pareto_k")
  throw_pareto_warnings(pareto_k)

  psis_object(
    unnormalized_log_weights = log_weights,
    pareto_k = pareto_k,
    tail_len = tail_len,
    r_eff = r_eff
  )
}

#' Extract named components from each list in the list of lists obtained by
#' parallelizing `do_psis_i()`
#'
#' @noRd
#' @param x List of lists.
#' @param item String naming the component or attribute to pull out of each list
#'   (or list-like object).
#' @param fun,fun.val passed to `vapply()`'s `FUN` and `FUN.VALUE` arguments.
#' @return Numeric vector or matrix.
#'
psis_apply <- function(x, item, fun = c("[[", "attr"), fun_val = numeric(1)) {
  if (!is.list(x)) stop("Internal error ('x' must be a list for psis_apply)")
  vapply(x, FUN = match.arg(fun), FUN.VALUE = fun_val, item)
}

#' PSIS on a single vector
#'
#' @noRd
#' @param log_ratios_i A vector of log importance ratios (for `loo()`, negative
#'   log likelihoods).
#' @param tail_len_i An integer tail length.
#'
#' @details
#' * The maximum of the log ratios is subtracted from each of them
#' * If there are enough tail samples then the tail is smoothed with PSIS
#' * The log weights (or log ratios if no smoothing) larger than zero are set to 0
#'
#' @return A named list containing:
#' * `lw`: vector of unnormalized log weights
#' * `pareto_k`: scalar Pareto k estimate.
#'
do_psis_i <- function(log_ratios_i, tail_len_i) {
  S <- length(log_ratios_i)
  lw_i <- log_ratios_i - max(log_ratios_i)
  khat <- Inf

  if (enough_tail_samples(tail_len_i)) {
    ord <- sort.int(lw_i, index.return = TRUE)
    tail_ids <- seq(S - tail_len_i + 1, S)
    lw_tail <- ord$x[tail_ids]
    if (abs(max(lw_tail) - min(lw_tail)) < sqrt(.Machine$double.eps)) {
      warning(
        "Can't fit generalized Pareto distribution ",
        "because all tail values are the same.",
        call. = FALSE
      )
    } else {
      cutoff <- ord$x[min(tail_ids) - 1] # largest value smaller than tail values
      smoothed <- psis_smooth_tail(lw_tail, cutoff)
      khat <- smoothed$k
      lw_i[ord$ix[tail_ids]] <- smoothed$tail
    }
  }

  # truncate at max of raw wts (i.e., 0 since max has been subtracted)
  lw_i[lw_i > 0] <- 0

  list(log_weights = lw_i, pareto_k = khat)
}


#' PSIS tail smoothing for a single vector
#'
#' @noRd
#' @param x Vector of tail elements already sorted in ascending order.
#' @return A named list containing:
#' * `tail`: vector same size as `x` containing the logs of the
#'   order statistics of the generalized pareto distribution.
#' * `k`: scalar shape parameter estimate.
#'
psis_smooth_tail <- function(x, cutoff) {
  len <- length(x)
  exp_cutoff <- exp(cutoff)

  # save time not sorting since x already sorted
  fit <- gpdfit(exp(x) - exp_cutoff, sort_x = FALSE)
  k <- fit$k
  sigma <- fit$sigma
  if (is.finite(k)) {
      p <- (seq_len(len) - 0.5) / len
      qq <- qgpd(p, k, sigma) + exp_cutoff
      tail <- log(qq)
  } else {
      tail <- x
  }
  list(tail = tail, k = k)
}


#' Calculate tail lengths to use for fitting the GPD
#'
#' The number of weights (i.e., tail length) used to fit the generalized Pareto
#' distribution is now decreasing with the number of posterior draws S, and is
#' also adjusted based on the relative MCMC neff for `exp(log_lik)`. This will
#' answer the questions about the asymptotic properties, works better for thick
#' tailed proposal distributions, and is adjusted based on dependent Markov chain
#' samples. Specifically, the tail length is now `3*sqrt(S)/r_eff` but capped at
#' 20% of the total number of weights.
#'
#' @noRd
#' @param r_eff A N-vector of relative MCMC effective sample sizes of `exp(log-lik matrix)`.
#' @param S The (integer) size of posterior sample.
#' @return An N-vector of tail lengths.
#'
n_pareto <- function(r_eff, S) {
  ceiling(pmin(0.2 * S, 3 * sqrt(S / r_eff)))
}

#' Check for enough tail samples to fit GPD
#'
#' @noRd
#' @param tail_len Integer tail length.
#' @param min_len The minimum allowed tail length.
#' @return `TRUE` or `FALSE`
#'
enough_tail_samples <- function(tail_len, min_len = 5) {
  tail_len >= min_len
}


#' Throw warnings about pareto k estimates
#'
#' @noRd
#' @param k A vector of pareto k estimates.
#' @param high The value at which to warn about slighly high estimates.
#' @param too_high The value at which to warn about very high estimates.
#' @return Nothing, just possibly throws warnings.
#'
throw_pareto_warnings <- function(k, high = 0.5, too_high = 0.7) {
  if (any(k > too_high)) {
    .warn("Some Pareto k diagnostic values are too high. ", .k_help())
  } else if (any(k > high)) {
    .warn("Some Pareto k diagnostic values are slightly high. ", .k_help())
  }
}

#' Warn if not enough tail samples to fit GPD
#'
#' @noRd
#' @param tail_lengths Vector of tail lengths.
#' @return `tail_lengths`, invisibly.
#'
throw_tail_length_warnings <- function(tail_lengths) {
  tail_len_bad <- !sapply(tail_lengths, enough_tail_samples)
  if (any(tail_len_bad)) {
    if (length(tail_lengths) == 1) {
      warning(
        "Not enough tail samples to fit the generalized Pareto distribution.",
        call. = FALSE, immediate. = TRUE
      )
    } else {
      bad <- which(tail_len_bad)
      Nbad <- length(bad)
      warning(
        "Not enough tail samples to fit the generalized Pareto distribution ",
        "in some or all columns of matrix of log importance ratios. ",
        "Skipping the following columns: ",
        paste(if (Nbad <= 10) bad else bad[1:10], collapse = ", "),
        if (Nbad > 10) paste0(", ... [", Nbad - 10, " more not printed].\n") else "\n",
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  invisible(tail_lengths)
}

#' Prepare `r_eff` to pass to `psis()` and throw warnings/errors if necessary
#'
#' @noRd
#' @param r_eff User's `r_eff` argument.
#' @param len The length `r_eff` should have if not `NULL` or `NA`.
#' @return
#' * If `r_eff` has length `len` then `r_eff` is returned.
#' * If `r_eff` is `NULL` then a warning is thrown and `rep(1, len)` is returned.
#' * If `r_eff` is `NA` then the warning is skipped and
#'   `rep(1, len)` is returned.
#' * If `r_eff` has length `len` but has `NA`s then an error is thrown.
#'
prepare_psis_r_eff <- function(r_eff, len) {
  if (isTRUE(is.null(r_eff) || all(is.na(r_eff)))) {
    if (!called_from_loo() && is.null(r_eff)) {
      throw_psis_r_eff_warning()
    }
    r_eff <- rep(1, len)
  } else if (length(r_eff) != len) {
    stop("'r_eff' must have one value per observation.", call. = FALSE)
  } else if (anyNA(r_eff)) {
    stop("Can't mix NA and not NA values in 'r_eff'.", call. = FALSE)
  }
  return(r_eff)
}

#' Check if `psis()` was called from one of the loo methods
#'
#' @noRd
#' @return `TRUE` if the `loo()` array, matrix, or function method is found in
#'   the active call list, `FALSE` otherwise.
#'
called_from_loo <- function() {
  calls <- sys.calls()
  txt <- unlist(lapply(calls, deparse))
  patts <- "loo.array\\(|loo.matrix\\(|loo.function\\("
  check <- sapply(txt, function(x) grepl(patts, x))
  isTRUE(any(check))
}

#' Warning message about missing `r_eff` argument
#' @noRd
throw_psis_r_eff_warning <- function() {
  warning(
    "Relative effective sample sizes ('r_eff' argument) not specified. ",
    "PSIS n_eff will not be adjusted based on MCMC n_eff.",
    call. = FALSE
  )
}

