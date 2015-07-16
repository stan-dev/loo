# Convenience functions ---------------------------------------------------

#' Named lists
#'
#' Create a named list using specified names or, if names are omitted, using the
#' names of the objects in the list. The code \code{list(a = a, b = b)} becomes
#' \code{nlist(a,b)} and \code{list(a = a, b = 2)} becomes \code{nlist(a, b =
#' 2)}, etc.
#'
#' @keywords internal
#' @export
#' @param ... Objects to include in the list.
#' @return A named list.
#'
#' @seealso \code{\link[base]{list}}
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
  out <- list(...)
  onms <- names(out)
  no_names <- is.null(onms)
  has_name <- if (no_names)
    FALSE else nzchar(onms)
  if (all(has_name))
    return(out)
  mc <- match.call(expand.dots = TRUE)
  nms <- as.character(mc)[-1L]
  if (no_names)
    names(out) <- nms
  else
    names(out)[!has_name] <- nms[!has_name]
  out
}

is.loo <- function(x) {
  inherits(x, "loo")
}

unlist_lapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...), use.names = FALSE)
}

cbind_list <- function(x) {
  do.call(cbind, x)
}
