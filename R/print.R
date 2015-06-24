#' Print methods
#' @param x a list of class \code{'loo'} (as returned by
#'   \code{\link{loo_and_waic}}) or of class \code{'compare.loo'} (as returned
#'   by \code{\link{loo_and_waic_diff}}).
#' @param ... Other arguments. Currently only \code{digits} is supported.
#' @export
print.loo <- function(x, ...) {
  dots <- list(...)
  digits <- 2
  if (length(dots) != 0 && "digits" %in% names(dots)) {
    digits <- dots$digits
  }
  dims <- attr(x, "log_lik_dim")
  L <- length(x)
  z <- x[-c(L - 1, L)]
  uz <- unlist(z)
  print_ord <- c(1, 3, 2, 4, 5, 6)
  out <- cbind(total = uz[print_ord], se = uz[print_ord + length(print_ord)])
  printCoefmat(out, digits = digits)
  cat("-----\n")
  cat(paste("Computed from", dims[1], "by", dims[2], "log-likelihood matrix"))
  invisible(x)
}


#' @rdname print.loo
#' @export
print.compare.loo <- function(x, ...) {
  dots <- list(...)
  digits <- 2
  if (length(dots) != 0 && "digits" %in% names(dots)) {
    digits <- dots$digits
  }
  ux <- unlist(x)
  print(round(ux, digits))
  invisible(x)
}
