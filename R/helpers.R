#' @importFrom matrixStats colLogSumExps
logColMeansExp <- function(x) {
  # should be more stable than log(colMeans(exp(x)))
  S <- nrow(x)
  colLogSumExps(x) - log(S)
}

# Inverse-CDF of generalized Pareto distribution
# (formula from Wikipedia)
qgpd <- function(p, xi = 1, mu = 0, beta = 1, lower.tail = TRUE) {
  if (!lower.tail)
    p <- 1 - p
  mu + beta * ((1 - p)^(-xi) - 1) / xi
}

# named lists
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names)
    FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }
  out
}

unlist_lapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...), use.names = FALSE)
}

seq_min_half <- function(L) {
  seq_len(L) - 0.5
}

vapply_seq <- function(L, FUN, ...) {
  vapply(1:L, FUN, FUN.VALUE = 0, ...)
}

cbind_list <- function(x) {
  stopifnot(is.list(x))
  do.call(cbind, x)
}

lx <- function(a, x) {
  k <- mean.default(log(1 - a * x))
  log(-a / k) - k - 1
}
