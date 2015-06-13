#' Print method for \code{'loo'} objects
#' @export
#' @param x a \code{'loo'} object returned by \code{\link{loo_and_waic}}.
#' @param ... other arguments. Currently only \code{digits} is supported.
#' @return Invisibly returns \code{x}.
#'
print.loo <- function(x, ...) {
  dots <- list(...)
  digits <- 2
  if (length(dots) != 0 && "digits" %in% names(dots)) {
    digts <- dots$digits
  }
  lx <- length(x)
  dims <- x[[lx]]
  z <- x[-c(lx-2, lx-1, lx)]
  uz <- unlist(z)
  ord <- c(1,3,2,4,5,6)
  out <- cbind(total = uz[ord], se = uz[ord + length(ord)])
  printCoefmat(out, digits = digits)
  cat("-----\n")
  cat(paste("Computed from", dims[1], "by", dims[2], "log_lik matrix"))
  invisible(x)
}

#' @export
#' @rdname print.loo
#'
print.compare.loo <- function(x, ...) {
  dots <- list(...)
  digits <- 2
  if (length(dots) != 0 && "digits" %in% names(dots)) {
    digts <- dots$digits
  }
  ux <- unlist(x)
  print(round(ux, digits))
  invisible(x)
}
