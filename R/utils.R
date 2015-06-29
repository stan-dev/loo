# Convenience functions ---------------------------------------------------

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

# cbind together the elements of a list
cbind_list <- function(x) {
  do.call(cbind, x)
}

# first n odd numbers
nodds <- function(n) {
  seq(1, by = 2, len = n)
}

# first n even numbers
nevens <- function(n) {
  seq(2, by = 2, len = n)
}
