# use matrixStats::colVars instead (much faster b/c written in C++)
# colVars <- function(X) {
#   # column variances
#   N <- dim(X)[[1]]
#   C <- dim(X)[[2]]
#   Xbar <- matrix(.colMeans(X, N, C), N, C, byrow = TRUE)
#   var <-  .colMeans((X - Xbar)^2, N, C) * N / (N - 1)
#   var
# }

nlist <- function(...) {
  # named lists
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name)) return(out)
  nms <- as.character(m)[-1]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }
  out
}

unlist_lapply <- function(x, fun) {
  unlist(lapply(x, fun), use.names = FALSE)
}

seq_min_half <- function(n) {
  seq_len(n) - 0.5
}

vapply1m <- function(m, fun) {
  vapply(1:m, fun, FUN.VALUE = 0)
}

lx <- function(a, x) {
  k <- mean.default(log(1 - a*x))
  log(-a/k) - k - 1
}

qgpd <- function(p, xi=1, mu=0, beta=1, lower.tail=TRUE){
  # Generalized Pareto inverse-cdf (formula from Wikipedia)
  if (!lower.tail) p <- 1-p
  mu + beta * ((1-p)^(-xi) - 1) / xi
}

sumlogs <- function(x, dimen=1) {
  # log_sum_exp
  x_max <- max(x)
  if (is.null(dim(x))) return(x_max + log(sum(exp(x-x_max)))) 
  if (dimen == 1) return(x_max + log(colSums(exp(x-x_max))))
  x_max + log(rowSums(exp(x-x_max)))
}

escape_check <- function(x, escape_if_greater_than = 700) {
  # find the difference between the largest and 2nd largest values in a matrix
  L <- length(x)
  x_sorted <- sort(x, method = "quick")[c(1:2, (L-1):L)]
  diffs <- vapply(c(1,3), function(i) abs(diff(x_sorted[i:(i+1)])), 0)
  if (any(diffs > escape_if_greater_than)) {
    message(paste0("Failed. Difference between largest and second largest", 
                   "log-weights is more than ", escape_if_greater_than,"."))
    return(invisible(NULL))
  }
}