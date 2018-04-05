#' Pareto smoothed importance sampling (PSIS)
#'
#' Implementation of Pareto smoothed importance sampling (PSIS), a method for
#' stabilizing importance ratios. The version of PSIS implemented here
#' corresponds to the algorithm presented in Vehtari, Gelman and Gabry (2017b).
#' For PSIS diagnostics see the \link{pareto-k-diagnostic} page.
#'
#' @export
#' @param log_ratios An array, matrix, or vector of importance ratios on the log
#'   scale (for PSIS-LOO these are \emph{negative} log-likelihood values). See the
#'   \strong{Methods (by class)} section below for a detailed description of how
#'   to specify the inputs for each method.
#' @param ... Arguments passed on to the various methods.
#' @template cores
#' @param r_eff Vector of relative effective sample size estimates containing
#'   one element per observation. The values provided should be the relative
#'   effective sample sizes of \code{1/exp(log_ratios)} (i.e., \code{1/ratios}).
#'   This is related to the relative efficiency of estimating the normalizing
#'   term in self-normalizing importance sampling. If \code{r_eff} is not
#'   provided then the reported PSIS effective sample sizes and Monte Carlo
#'   error estimates will be over-optimistic. See the \code{\link{relative_eff}}
#'   helper function for computing \code{r_eff}.
#'
#' @return The \code{psis} methods return an object of class \code{"psis"},
#'   which is a named list with the following components:
#'
#' \describe{
#'   \item{\code{log_weights}}{
#'     Vector or matrix of smoothed (and truncated) but \emph{unnormalized} log
#'     weights. To get normalized weights use the \code{weights} method provided
#'     for objects of class \code{"psis"}.
#'   }
#'  \item{\code{diagnostics}}{
#'    A named list containing two vectors:
#'    \itemize{
#'     \item \code{pareto_k}: Estimates of the shape parameter \eqn{k} of the
#'     generalized Pareto distribution. See the \code{\link{pareto-k-diagnostic}}
#'     page for details.
#'     \item \code{n_eff}: PSIS effective sample size estimates.
#'   }
#'  }
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
#'   \item{\code{r_eff}}{
#'     If specified, the user's \code{r_eff} argument.
#'   }
#'   \item{\code{dims}}{
#'     Integer vector of length 2 containing \code{S} (posterior sample size)
#'     and \code{N} (number of observations).
#'   }
#' }
#'
#' @seealso
#' \itemize{
#' \item \code{\link{loo}} for approximate LOO-CV using PSIS.
#' \item \code{\link{pareto-k-diagnostic}} for PSIS diagnostics.
#' }
#'
#' @template loo-and-psis-references
#'
#' @examples
#' log_ratios <- -1 * example_loglik_array()
#' r_eff <- relative_eff(exp(log_ratios))
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
#'
psis <- function(log_ratios, ...) UseMethod("psis")

#' @export
#' @templateVar fn psis
#' @template array
#'
psis.array <-
  function(log_ratios, ..., r_eff = NULL, cores = getOption("mc.cores", 1)) {
    cores <- loo_cores(cores)
    stopifnot(
      length(dim(log_ratios)) == 3,
      is.null(r_eff) || length(r_eff) == dim(log_ratios)[3]
    )
    log_ratios <- validate_ll(log_ratios)
    log_ratios <- llarray_to_matrix(log_ratios)
    if (is.null(r_eff)) {
      if (!called_from_loo()) throw_psis_r_eff_warning()
      r_eff <- rep(1, ncol(log_ratios))
    }

    do_psis(log_ratios, r_eff = r_eff, cores = cores)
  }

#' @export
#' @templateVar fn psis
#' @template matrix
#'
psis.matrix <-
  function(log_ratios, ..., r_eff = NULL, cores = getOption("mc.cores", 1)) {
    cores <- loo_cores(cores)
    stopifnot(is.null(r_eff) || length(r_eff) == ncol(log_ratios))
    log_ratios <- validate_ll(log_ratios)
    if (is.null(r_eff)) {
      if (!called_from_loo()) throw_psis_r_eff_warning()
      r_eff <- rep(1, ncol(log_ratios))
    }
    do_psis(log_ratios, r_eff = r_eff, cores = cores)
  }

#' @export
#' @templateVar fn psis
#' @template vector
#'
psis.default <-
  function(log_ratios, ..., r_eff = NULL) {
    stopifnot(is.null(dim(log_ratios)) || length(dim(log_ratios)) == 1,
              is.null(r_eff) || length(r_eff) == 1)
    dim(log_ratios) <- c(length(log_ratios), 1)
    if (is.null(r_eff)) {
      if (!called_from_loo()) throw_psis_r_eff_warning()
      r_eff <- 1
    }
    psis.matrix(log_ratios,
                r_eff = r_eff,
                cores = 1)
  }

#' @rdname psis
#' @export
#' @export weights.psis
#' @method weights psis
#' @param object For the \code{weights} method, an object
#'   returned by \code{psis} (a list with class \code{"psis"}).
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


#' @export
dim.psis <- function(x) {
  attr(x, "dims")
}


# internal ----------------------------------------------------------------

# Do PSIS given matrix of log weights
#
# @param lr Matrix of log ratios (-loglik)
# @param r_eff Vector of relative effective sample sizes
# @param cores User's integer cores argument
# @return A list with class "psis" and structure described in the main doc at
#   the top of this file.
#
do_psis <- function(log_ratios, r_eff, cores) {
  stopifnot(length(r_eff) == ncol(log_ratios),
            cores == as.integer(cores))
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
        FUN = function(i)
          do_psis_i(log_ratios[, i], tail_len[i]),
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
            do_psis_i(log_ratios[, i], tail_len[i])
        )
    }
  }

  pareto_k <- psis_apply(lw_list, "pareto_k")
  throw_pareto_warnings(pareto_k)

  log_weights <- psis_apply(lw_list, "log_weights", fun_val = numeric(S))

  # truncate at the max of raw wts
  max_wts <- apply(log_ratios, 2, max)
  log_weights <- sapply(1:N, function(n) {
    truncate_lw(log_weights[, n], upper = max_wts[n])
  })

  psis_object(
    unnormalized_log_weights = log_weights,
    pareto_k = pareto_k,
    tail_len = tail_len,
    r_eff = r_eff
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
#' @param r_eff Vector of relative MCMC n_eff for exp(log lik)
#' @return A list of class "psis" with structure described in the main doc at the
#'   top of this file.
#'
psis_object <-
  function(unnormalized_log_weights,
           pareto_k,
           tail_len,
           r_eff) {
    stopifnot(is.matrix(unnormalized_log_weights))

    norm_const_log <- colLogSumExps(unnormalized_log_weights)
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
    w <- weights(out, normalize = TRUE, log = FALSE)
    out$diagnostics[["n_eff"]] <- psis_n_eff(w, r_eff)
    return(out)
  }


#' PSIS (without truncation and normalization) on a single vector
#'
#' @noRd
#' @param log_ratios_i A vector of log ratios (for loo negative log likelihoods).
#' @param tail_len_i An integer tail length.
#' @return list containing lw (vector of untruncated and unnormalized log
#'   weights), and pareto_k ('khat' estimate).
#'
do_psis_i <- function(log_ratios_i, tail_len_i) {
  S <- length(log_ratios_i)
  lw_i <- log_ratios_i - max(log_ratios_i)
  khat <- Inf

  if (enough_tail_samples(tail_len_i)) {
    ord <- sort.int(lw_i, index.return = TRUE)
    tail_ids <- seq(S - tail_len_i + 1, S)
    lw_tail <- ord$x[tail_ids]
    if (all(lw_tail == lw_tail[1])) {
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
  list(log_weights = lw_i, pareto_k = khat)
}


#' PSIS tail smoothing for a single vector
#'
#' @noRd
#' @param x Vector of tail elements already sorted in ascending order.
#' @return List with components 'tail' (vector same size as x) and 'k' (scalar
#'   shape parameter estimate). The 'tail' component contains the logs of the
#'   order statistics of the GPD.
#'
psis_smooth_tail <- function(x, cutoff) {
  len <- length(x)
  exp_cutoff <- exp(cutoff)

  # save time not sorting since x already sorted
  fit <- gpdfit(exp(x) - exp_cutoff, sort_x = FALSE)
  k <- fit$k
  sigma <- fit$sigma

  p <- (seq_len(len) - 0.5) / len
  qq <- qgpd(p, k, sigma) + exp_cutoff
  list(tail = log(qq), k = k)
}


#' Calculate tail lengths to use for fitting the GPD
#'
#' The number of weights (i.e., tail length) used to fit the generalized Pareto
#' distribution is now decreasing with the number of posterior draws S, and is
#' also adjusted based on the relative MCMC neff for exp(log_lik). This will
#' answer the questions about the asymptotic properties, works better for thick
#' tailed proposal distributions, and is adjusted based on dependent Markov chain
#' samples. Specifically, the tail length is now 3*sqrt(S)/r_eff but capped at
#' 20% of the total number of weights.
#'
#' @noRd
#' @param r_eff A N-vector of relative MCMC effective sample sizes of
#'   exp(log-lik matrix).
#' @param S The (integer) size of posterior sample.
#' @return An N-vector of tail lengths.
#'
n_pareto <- function(r_eff, S) {
  n <- pmin(0.2 * S, 3 * sqrt(S / r_eff))
  ceiling(n)
}

#' Check for enough tail samples to fit GPD
#'
#' @noRd
#' @param tail_len Integer tail length.
#' @param min_len The minimum allowed tail length.
#' @return TRUE or FALSE
#'
enough_tail_samples <- function(tail_len, min_len = 5) {
  tail_len >= min_len
}


#' Truncate a vector of (log) weights
#'
#' @noRd
#' @param x A vector of log weights.
#' @param upper Upper bound.
#' @return x, but modified so all values above \code{upper} are set equal to \code{upper}.
#'
truncate_lw <- function(x, upper) {
  x[x > upper] <- upper
  return(x)
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
    .warn("Some Pareto k diagnostic values are too high. ",
          .k_help())
  } else if (any(k > high)) {
    .warn("Some Pareto k diagnostic values are slightly high. ",
          .k_help())
  }
}

#' Warn if not enough tail samples to fit GPD
#'
#' @noRd
#' @param tail_lengths Vector of tail lengths.
#' @return tail_lengths, invisibly.
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
        call. = FALSE, immediate. = TRUE
      )
    }
  }
  invisible(tail_lengths)
}


#' Warning message if r_eff not specified
#' @noRd
#'
throw_psis_r_eff_warning <- function() {
  warning("Relative effective sample sizes ('r_eff' argument) not specified. ",
          "PSIS n_eff will not adjusted based on MCMC n_eff.", call. = FALSE)
}

#' check if psis was called from one of the loo methods
#' @noRd
called_from_loo <- function() {
  calls <- sys.calls()
  txt <- unlist(lapply(calls, deparse))
  patts <- "loo.array\\(|loo.matrix\\(|loo.function\\("
  check <- sapply(txt, function(x) grepl(patts, x))
  isTRUE(any(check))
}

