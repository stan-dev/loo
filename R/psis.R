#' Pareto smoothed importance sampling (PSIS)
#'
#' Implementation of Pareto smoothed importance sampling, a method for
#' stabilizing importance weights. For full details about the algorithm see
#' Vehtari, Gelman and Gabry (2016, 2017). For diagnostics see the
#' \link{pareto-k-diagnostic} page.
#'
#' @export
#' @param x A log-likelihood array, matrix, or function. See the \strong{Methods
#'   (by class)} section below for a detailed description.
#'
#' @param args Only required if \code{x} is a function. A list containing
#'   the data required to specify the arguments to the function. See the
#'   \strong{Methods (by class)} section below for how \code{args} should be
#'   specified.
#'
#' @return \code{psis} returns a named list of class "psis" with components
#'   \code{lw_smooth} (smoothed but \emph{unnormalized} log weights) and
#'   \code{pareto_k} (estimated generalized Pareto
#'   \link[=pareto-k-diagnostic]{shape parameter(s) k}). To get normalized log
#'   weights use the \code{weights} method for objects of class "psis".
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
psis <- function(x, ...) UseMethod("psis")

#' @export
#' @describeIn psis
#' An iterations by chains by observations array.
#'
psis.array <- function(x, ..., wtrunc = 3/4, cores = getOption("loo.cores", 1)) {
  stopifnot(length(dim(x)) == 3, !anyNA(x))
  rel_neff <- relative_neff(exp(x))

  xdim <- dim(x)
  dim(x) <- c(xdim[1] * xdim[2], xdim[3])
  x <- unname(x)
  n_tail <- n_pareto(rel_neff, S = nrow(x))
  do_psis(
    lw = -x,
    tail_len = n_tail,
    wtrunc = wtrunc,
    cores = cores
  )
}

#' @export
#' @describeIn psis
#' An iterations by observations matrix.
#'
psis.matrix <- function(x, chain_id = NULL, ..., wtrunc = 3/4, cores = getOption("loo.cores", 1)) {
  stopifnot(!anyNA(x))
  rel_neff <- relative_neff(exp(x), chain_id)
  n_tail <- n_pareto(rel_neff, S = nrow(x))
  do_psis(
    lw = -x,
    tail_len = n_tail,
    wtrunc = wtrunc,
    cores = cores
  )
}

#' @export
#' @describeIn psis
#' A vector.
psis.vector <- function(x, chain_id = NULL, ..., wtrunc = 3/4) {
  stopifnot(!anyNA(x))
  dim(x) <- c(length(x), 1)
  psis.matrix(x,
              chain_id = chain_id,
              wtrunc = wtrunc,
              cores = 1)
}

#' @export
#' @templateVar fn psis
#' @template function
#'
psis.function <- function(x, args, ..., wtrunc = 3/4) {
  if (missing(args))
    stop("'args' must be specified.")
}


# internal ----------------------------------------------------------------

# Do PSIS given matrix of log weights
#
# @param lw matrix of log weights (-loglik)
# @param vector of tail lengths
# @param user's scalar wtrunc parameter
#
do_psis <- function(lw, tail_len, wtrunc, cores) {
  stopifnot(is.matrix(lw), length(tail_len) == ncol(lw))
  S <- nrow(lw)
  N <- ncol(lw)

  if (cores == 1) {
    lw_list <- lapply(seq_len(N), function(i)
      .psis_i(lw[, i], tail_len[i]))
  } else {
    if (.Platform$OS.type != "windows") {
      lw_list <- mclapply(
        X = seq_len(N),
        FUN = function(i) .psis_i(lw[, i], tail_len[i]),
        mc.cores = cores
      )
    } else { # nocov start
      cl <- makePSOCKcluster(cores)
      on.exit(stopCluster(cl))
      lw_list <-
        parLapply(
          cl = cl,
          X = seq_len(N),
          fun = function(i) .psis_i(lw[, i], tail_len[i])
        )
    } # nocov end
  }

  lw_smooth <- vapply(lw_list, "[[", "lw", FUN.VALUE = numeric(S))
  if (wtrunc > 0) {
    lw_smooth <- apply(lw_smooth, 2, "truncate_lw",
                       logS = log(S), wtrunc = wtrunc)
  }

  pareto_k <- vapply(lw_list, "[[", "pareto_k", FUN.VALUE = numeric(1))
  throw_psis_warnings(pareto_k)

  structure(
    nlist(lw_smooth, pareto_k),
    const = matrixStats::colLogSumExps(lw_smooth),
    tail_len = tail_len,
    wtrunc = wtrunc,
    class = c("psis", "list")
  )
}

# PSIS (without truncation and normalization) on a single vector
#
# @param lw_i A vector of log weights.
# @param tail_len_i An integer tail length.
# @return list containing lw (vector of untruncated and unnormalized log
#   weights), and pareto_k ('khat' estimate).
#
.psis_i <- function(lw_i, tail_len_i) {
  S <- length(lw_i)
  lw_i <- lw_i - max(lw_i)

  if (!enough_tail_samples(tail_len_i)) {
    warning(
      "Too few tail samples to fit generalized Pareto distribution.\n",
      "Weights will be truncated but not smoothed.",
      call. = FALSE
    )
    lw_new <- lw_i
    k <- Inf
  } else {
    tail_ids <- seq(S - tail_len_i + 1, S)
    ord <- sort.int(lw_i, index.return = TRUE)
    lw_tail <- ord$x[tail_ids]
    if (all(lw_tail == lw_tail[1])) {
      warning(
        "All tail values are the same.\n",
        "Weights are truncated but not smoothed.",
        call. = FALSE
      )
      lw_new <- lw_i
      k <- Inf
    } else {
      smoothed <- psis_smooth_tail(lw_tail)
      ord$x[tail_ids] <- smoothed$tail
      k <- smoothed$k
      lw_new <- ord$x[ord$ix]
    }
  }

  structure(
    list(lw = lw_new, pareto_k = k),
    class = c("psis_i", "list")
  )
}

# calculate tail length for PSIS
n_pareto <- function(rel_neff, S) {
  ceiling(pmin(0.2 * S, 3 * sqrt(S) / rel_neff))
}

# @param x Vector of sorted tail elements
psis_smooth_tail <- function(x) {
  tail_len <- length(x)
  exp_cutoff <- exp(x[1])
  fit <- gpdfit(exp(x) - exp_cutoff)
  prb <- (seq_len(tail_len) - 0.5) / tail_len
  qq <- qgpd(p = prb, xi = fit$k, sigma = fit$sigma) + exp_cutoff
  list(tail = log(qq), k = fit$k)
}

# @param x vector of log weights
# @param logS precomputed log(length(x))
# @param wtrunc user's wtrunc parameter
truncate_lw <- function(x, logS, wtrunc) {
  lwtrunc <- wtrunc * logS - logS + logSumExp(x)
  x[x > lwtrunc] <- lwtrunc
  return(x)
}

# inverse-CDF of generalized Pareto distribution (formula from Wikipedia)
qgpd <-
  function(p,
           xi = 1,
           mu = 0,
           sigma = 1,
           lower.tail = TRUE) {
    if (is.nan(sigma) || sigma <= 0)
      return(rep(NaN, length(p)))
    if (!lower.tail)
      p <- 1 - p

    mu + sigma * ((1 - p) ^ (-xi) - 1) / xi
  }

enough_tail_samples <- function(tail_len, min_len = 5) {
  tail_len >= min_len
}


# warnings about pareto k values ------------------------------------------
throw_psis_warnings <- function(k) {
  if (any(k > 0.7)) {
    .warn("Some Pareto k diagnostic values are too high. ", .k_help())
  } else if (any(k > 0.5)) {
    .warn("Some Pareto k diagnostic values are slightly high. ", .k_help())
  }
}

#' @rdname psis
#' @export
#' @method weights psis
#' @param object For the \code{weights} method, an object returned by
#'   \code{psis}, which is a list with class "psis".
#' @param log For the \code{weights} method, should the weights be returned on
#'   the log scale? Defaults to \code{TRUE}.
#' @param normalize For the \code{weights} method, should the weights be
#'   normalized? Defaults to \code{TRUE}.
#'
#' @return The \code{weights} method returns a matrix, which by default is a
#'   matrix of normalized log weights.
#'
weights.psis <-
  function(object,
           ...,
           log = TRUE,
           normalize = TRUE) {
    out <- object[["lw_smooth"]]
    const <- attr(object, "const") # precomputed colLogSumExp(lw)
    if (normalize)
      out <- sweep(out, 2, const)
    if (!log)
      out <- exp(out)

    return(out)
  }
