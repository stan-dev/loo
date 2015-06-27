#' @importFrom matrixStats colVars
.pointwise_waic <- function(log_lik) {
  lpd <- logColMeansExp(log_lik)
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- .ic(elpd_waic)
  nlist(elpd_waic, p_waic, waic)
}
.pointwise_loo <- function(log_lik, vgis) {
  # vgis is output from vgisloo()
  lpd <- logColMeansExp(log_lik)
  elpd_loo <- vgis$loos
  p_loo <- lpd - elpd_loo
  looic <- .ic(elpd_loo)
  nlist(elpd_loo, p_loo, looic)
}
.totals <- function(pointwise) {
  N <- length(pointwise[[1L]])
  total  <- unlist_lapply(pointwise, sum)
  se <- sqrt(N * unlist_lapply(pointwise, var))
  as.list(c(total, se))
}
.ic <- function(elpd) {
  stopifnot(is.vector(elpd), is.numeric(elpd))
  -2 * elpd
}

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

lx <- function(a, x) {
  k <- mean.default(log1p(-a * x))
  log(-a / k) - k - 1
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

is.loo <- function(x) {
  inherits(x, "loo")
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
  # cbind together the elements of a list
  do.call(cbind, x)
}

c_list <- function(x) {
  # concatenate elements of a list into a vector
  do.call(c, x)
}

# first n odd numbers
nodds <- function(n) {
  seq(1, by = 2, len = n)
}

# first n even numbers
nevens <- function(n) {
  seq(2, by = 2, len = n)
}
