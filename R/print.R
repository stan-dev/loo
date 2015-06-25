#' Print methods
#' @export
#'
#' @param x a list of class \code{'loo'} (as returned by
#'   \code{\link{loo_and_waic}}) or of class \code{'compare.loo'} (as returned
#'   by \code{\link{loo_and_waic_diff}}).
#' @param ... ignored.
#' @param digits number of significant digits to display.
#' @return returns \code{x} invisibly.
#' @seealso \code{\link{loo-package}}, \code{\link{loo_and_waic}},
#' \code{\link{loo_and_waic_diff}}
#'
print.loo <- function(x, ..., digits = 1) {
  dims <- attr(x, "log_lik_dim")
  z <- x[-grep("pointwise|pareto_k", names(x))]
  uz <- unlist(z)
  nms <- names(uz)
  ses <- grepl("se", nms)
  out <- cbind(Estimate = uz[!ses], StdError = uz[ses])
  out <- format(round(out, digits), nsmall = digits)
  print(out, quote = FALSE)
  cat("-----\n")
  cat(paste("Computed from", dims[1], "by", dims[2], "log-likelihood matrix"))
  invisible(x)
}


#' @rdname print.loo
#' @export
print.compare.loo <- function(x, ..., digits = 1) {
  ux <- unlist(x)
  out <- format(round(ux, digits), nsmall = digits)
  print(out, quote = FALSE)
  invisible(x)
}
