#' Print methods
#' @export
#'
#' @param x a list of class \code{'loo'} (as returned by \code{\link{loo}} or
#'   \code{\link{waic}}) or of class \code{'compare.loo'} (as returned by
#'   \code{\link{compare}}).
#' @param ... ignored.
#' @param digits number of significant digits to display.
#' @return Returns \code{x} invisibly.
#'
print.loo <- function(x, ..., digits = 1) {
  dims <- attr(x, "log_lik_dim")
  msg <- paste("Computed from", dims[1], "by", dims[2],
               "log-likelihood matrix\n\n")
  cat(msg)
  z <- x[-grep("pointwise|pareto_k", names(x))]
  uz <- unlist(z)
  nms <- names(uz)
  ses <- grepl("se", nms)
  out <- data.frame(Estimate = uz[!ses], SE = uz[ses])
  out <- format(round(out, digits), nsmall = digits)
  print(out, quote = FALSE)
  invisible(x)
}

#' @rdname print.loo
#' @export
print.compare.loo <- function(x, ..., digits = 1) {
  ux <- unlist(x)
  names(ux) <- c("elpd_diff", "SE")
  out <- format(round(ux, digits), nsmall = digits)
  print(out, quote = FALSE)
  invisible(x)
}
