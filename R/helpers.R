# waic and loo helpers ----------------------------------------------------

#' @importFrom matrixStats colLogSumExps
logColMeansExp <- function(x) {
  # should be more stable than log(colMeans(exp(x)))
  S <- nrow(x)
  colLogSumExps(x) - log(S)
}

logColMeansExp_ll <- function(fun, args) {
  # should be more stable than log(colMeans(exp(x)))
  logS <- log(args$S)
  clse <- vapply(seq_len(args$N), FUN = function(i) {
    logSumExp(fun(i = i, data = args$data[i,,drop=FALSE], draws = args$draws))
  }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
  clse - logS
}

colVars_ll <- function(fun, args) {
  vapply(seq_len(args$N), FUN = function(i) {
    x <- fun(i = i, data = args$data[i,,drop=FALSE], draws = args$draws)
    var(as.vector(x))
  }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
}

totals <- function(pointwise) {
  N <- length(pointwise[[1L]])
  total  <- unlist_lapply(pointwise, sum)
  se <- sqrt(N * unlist_lapply(pointwise, var))
  as.list(c(total, se))
}

#' @importFrom matrixStats colVars
pointwise_waic <- function(log_lik, llfun = NULL, llargs = NULL) {
  if (!missing(log_lik)) {
    lpd <- logColMeansExp(log_lik)
    p_waic <- colVars(log_lik)
  } else {
    if (is.null(llfun) || is.null(llargs))
      stop("Either 'log_lik' or 'llfun' and 'llargs' must be specified.",
           call. = FALSE)
    lpd <- logColMeansExp_ll(llfun, llargs)
    p_waic <- colVars_ll(llfun, llargs)
  }
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  pointwise <- nlist(elpd_waic, p_waic, waic)
  out <- totals(pointwise)
  nms <- names(pointwise)
  names(out) <- c(nms, paste0("se_", nms))
  out$pointwise <- cbind_list(pointwise)
  out
}
pointwise_loo <- function(psis, log_lik, llfun = NULL, llargs = NULL) {
  if (!missing(log_lik)) {
    lpd <- logColMeansExp(log_lik)
  } else {
    if (is.null(llfun) || is.null(llargs))
      stop("Either 'log_lik' or 'llfun' and 'llargs' must be specified.",
           call. = FALSE)
    lpd <- logColMeansExp_ll(llfun, llargs)
  }
  elpd_loo <- psis$loos
  p_loo <- lpd - elpd_loo
  looic <- -2 * elpd_loo
  pointwise <- nlist(elpd_loo, p_loo, looic)
  out <- totals(pointwise)
  nms <- names(pointwise)
  names(out) <- c(nms, paste0("se_", nms))
  out$pointwise <- cbind_list(pointwise)
  out$pareto_k <- psis$pareto_k
  out
}

# print and warning helpers -----------------------------------------------
.fr <- function(x, digits) format(round(x, digits), nsmall = digits)
.warn <- function(..., call. = FALSE) warning(..., call. = call.)
.k_help <- function() "See help('pareto-k-diagnostic') for details."
.k_cut <- function(k) {
  cut(
    k,
    breaks = c(-Inf, 0.5, 0.7, 1, Inf),
    labels = c("(-Inf, 0.5]", "(0.5, 0.7]", "(0.7, 1]", "(1, Inf)")
  )
}

pwaic_warnings <- function(p, digits = 1) {
  badp <- p > 0.4
  if (any(badp)) {
    count <- sum(badp)
    prop <- count / length(badp)
    .warn(paste0(count, " (", .fr(100 * prop, digits),
                 "%) p_waic estimates greater than 0.4."),
          "\nWe recommend trying loo() instead.")
  }
  invisible(NULL)
}


# convenience functions ---------------------------------------------------
is.loo <- function(x) {
  inherits(x, "loo")
}
unlist_lapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...), use.names = FALSE)
}
cbind_list <- function(x) {
  do.call(cbind, x)
}

#' Named lists
#'
#' Create a named list using specified names or, if names are omitted, using the
#' names of the objects in the list. The code \code{list(a = a, b = b)} becomes
#' \code{nlist(a,b)} and \code{list(a = a, b = 2)} becomes \code{nlist(a, b =
#' 2)}, etc.
#'
#' @export
#' @keywords internal
#' @param ... Objects to include in the list.
#' @return A named list.
#'
#' @seealso \code{\link[base]{list}}
#' @author Jonah Gabry
#' @examples
#'
#' # All variables already defined
#' a <- rnorm(100)
#' b <- mat.or.vec(10, 3)
#' nlist(a,b)
#'
#' # Define some variables in the call and take the rest from the environment
#' nlist(a, b, veggies = c("lettuce", "spinach"), fruits = c("banana", "papaya"))
#'
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }

  return(out)
}

# nocov start
# release reminders (for devtools)
release_questions <- function() {
  c(
    "Have you updated all references to the LOO paper?",
    "Have you updated inst/CITATION?",
    "Have you updated R code in vignette to match the code in the paper?"
  )
}
# nocov end
