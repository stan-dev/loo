#' @rdname loo_and_waic
#' @export
print.loo <- function(x, ...) {
  dots <- list(...)
  digits <- 2
  if (length(dots) != 0 && "digits" %in% names(dots)) {
    digits <- dots$digits
  }
  L <- length(x)
  dims <- x[[L]]
  z <- x[-c(L - 2, L - 1, L)]
  uz <- unlist(z)
  ord <- c(1, 3, 2, 4, 5, 6)
  out <- cbind(total = uz[ord], se = uz[ord + length(ord)])
  printCoefmat(out, digits = digits)
  cat("-----\n")
  cat(paste("Computed from", dims[1], "by", dims[2], "log_lik matrix"))
  invisible(x)
}


#' @rdname loo_and_waic_diff
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
