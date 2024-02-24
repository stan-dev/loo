#' Pareto smoothed importance sampling (PSIS)
#'
#' Implementation of Pareto smoothed importance sampling (PSIS), a method for
#' stabilizing importance ratios. The version of PSIS implemented here
#' corresponds to the algorithm presented in Vehtari, Simpson, Gelman, Yao,
#' and Gabry (2022).
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
#'   error estimates can be over-optimistic. If the posterior draws are (near)
#'   independent then `r_eff=1` can be used. `r_eff` has to be a scalar (same
#'   value is used for all observations) or a vector with length equal to the
#'   number of observations. The default value is 1. See the [relative_eff()]
#'   helper function for computing `r_eff`.
#'
#' @return The `psis()` methods return an object of class `"psis"`,
#'   which is a named list with the following components:
#'
#' \describe{
#'   \item{`log_weights`}{
#'     Vector or matrix of smoothed (and truncated) but *unnormalized* log
#'     weights. To get normalized weights use the
#'     [`weights()`][weights.importance_sampling] method provided for objects of
#'     class `"psis"`.
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
#'   \item{`method`}{
#'     Method used for importance sampling, here `psis`.
#'   }
#' }
#'
#' @seealso
#' * [loo()] for approximate LOO-CV using PSIS.
#' * [pareto-k-diagnostic] for PSIS diagnostics.
#' * The __loo__ package [vignettes](https://mc-stan.org/loo/articles/index.html)
#'   for demonstrations.
#' * The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html) on
#'   the __loo__ website for answers to frequently asked questions.
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
#'
psis <- function(log_ratios, ...) UseMethod("psis")

#' @export
#' @templateVar fn psis
#' @template array
#'
psis.array <-
  function(log_ratios, ...,
           r_eff = 1,
           cores = getOption("mc.cores", 1)) {
  importance_sampling.array(log_ratios = log_ratios, ...,
                            r_eff = r_eff,
                            cores = cores,
                            method = "psis")
  }


#' @export
#' @templateVar fn psis
#' @template matrix
#'
psis.matrix <-
  function(log_ratios,
           ...,
           r_eff = 1,
           cores = getOption("mc.cores", 1)) {
    importance_sampling.matrix(log_ratios,
                               ...,
                               r_eff = r_eff,
                               cores = cores,
                               method = "psis")
  }

#' @export
#' @templateVar fn psis
#' @template vector
#'
psis.default <-
  function(log_ratios, ..., r_eff = 1) {
    importance_sampling.default(log_ratios = log_ratios, ...,
                                r_eff = r_eff,
                                method = "psis")
  }


#' @rdname psis
#' @export
#' @param x For `is.psis()`, an object to check.
is.psis <- function(x) {
  inherits(x, "psis") && is.list(x)
}


# internal ----------------------------------------------------------------

#' @noRd
#' @seealso importance_sampling_object
psis_object <-
  function(unnormalized_log_weights,
           pareto_k,
           tail_len,
           r_eff) {
    importance_sampling_object(unnormalized_log_weights = unnormalized_log_weights,
                               pareto_k = pareto_k,
                               tail_len = tail_len,
                               r_eff = r_eff,
                               method = "psis")
  }


#' @noRd
#' @seealso do_importance_sampling
do_psis <- function(log_ratios, r_eff, cores, method){
  do_importance_sampling(log_ratios = log_ratios,
                         r_eff = r_eff,
                         cores = cores,
                         method = "psis")
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
#' @param ... Not used. Included to conform to API for differen IS methods.
#'
#' @details
#' * If there are enough tail samples then the tail is smoothed with PSIS
#' * The log weights (or log ratios if no smoothing) larger than the largest raw
#'   ratio are set to the largest raw ratio
#'
#' @return A named list containing:
#' * `lw`: vector of unnormalized log weights
#' * `pareto_k`: scalar Pareto k estimate.
#'
do_psis_i <- function(log_ratios_i, tail_len_i, ...) {
  S <- length(log_ratios_i)
  # shift log ratios for safer exponentation
  lw_i <- log_ratios_i - max(log_ratios_i)
  khat <- Inf

  if (enough_tail_samples(tail_len_i)) {
    ord <- sort.int(lw_i, index.return = TRUE)
    tail_ids <- seq(S - tail_len_i + 1, S)
    lw_tail <- ord$x[tail_ids]
    if (abs(max(lw_tail) - min(lw_tail)) < .Machine$double.eps/100) {
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
  # shift log weights back so that the smallest log weights remain unchanged
  lw_i <- lw_i + max(log_ratios_i)

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
#' @param r_eff A N-vector or scalar of relative MCMC effective sample sizes of
#'   `exp(log-lik matrix)`. The default value is 1.
#' @param S The (integer) size of posterior sample.
#' @return An N-vector of tail lengths.
#'
n_pareto <- function(r_eff, S) {
  if (isTRUE(is.null(r_eff) || all(is.na(r_eff)))) {
    r_eff <- 1
  }
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


#' Throw warnings about Pareto k estimates
#'
#' @noRd
#' @param k A vector of Pareto k estimates.
#' @param k_threshold The value at which to warn about high Pareto k estimates.
#' @return Nothing, just possibly throws warnings.
#'
throw_pareto_warnings <- function(k, k_threshold) {
  if (isTRUE(any(k > k_threshold))) {
    .warn("Some Pareto k diagnostic values are too high. ", .k_help())
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
#' * If `r_eff` is `NULL` then `rep(1, len)` is returned.
#' * If `r_eff` is `NA` then `rep(1, len)` is returned.
#' * If `r_eff` is a scalar then `rep(r_eff, len)` is returned.
#' * If `r_eff` is not a scalar but the length is not `len` then an error is thrown.
#' * If `r_eff` has length `len` but has `NA`s then an error is thrown.
#'
prepare_psis_r_eff <- function(r_eff, len) {
  if (isTRUE(is.null(r_eff) || all(is.na(r_eff)))) {
    r_eff <- rep(1, len)
  } else if (length(r_eff) == 1) {
    r_eff <- rep(r_eff, len)
  } else if (length(r_eff) != len) {
    stop("'r_eff' must have one value or one value per observation.", call. = FALSE)
  } else if (anyNA(r_eff)) {
    stop("Can't mix NA and not NA values in 'r_eff'.", call. = FALSE)
  }
  r_eff
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
    "PSIS ESS (n_eff) will not be adjusted based on MCMC ESS (n_eff).",
    call. = FALSE
  )
}

