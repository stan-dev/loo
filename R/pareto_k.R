#' Identify problematic observations
#'
#' Print a diagnostic table summarizing the estimated Pareto shape parameters or
#' find the indexes of observations for which the estimated Pareto shape
#' parameter \eqn{k} is larger than some \code{threshold} value.
#'
#' @name pareto-k-diagnostic
#' @param x For \code{pareto_k_ids} and \code{pareto_k_table}, an object created
#'   by \code{\link{loo}} or \code{\link{psislw}}. For the print method for
#'   class "pareto_k_table" \code{x} should be the object returned by
#'   \code{pareto_k_table}.
#' @details  See the PSIS-LOO section in \code{\link{loo-package}} for details
#'   about the interpretation of \eqn{k}.
#'
#' @return
#' \code{pareto_k_table} returns an object of class \code{"pareto_k_table"},
#' which is a matrix with columns \code{"Count"} and \code{"Proportion"} and has
#' its own print method.
#'
#' \code{pareto_k_ids} returns an integer vector indicating which observations
#' have Pareto \eqn{k} estimates above \code{threshold}.
#'
NULL

#' @rdname pareto-k-diagnostic
#' @export
#' @param threshold The threshold value for \eqn{k}.
pareto_k_ids <- function(x, threshold = 0.5) {
  k <- get_pareto_k(x)
  which(k > threshold)
}

#' @rdname pareto-k-diagnostic
#' @export
pareto_k_table <- function(x) {
  k <- get_pareto_k(x)
  kcut <- .k_cut(k)
  count <- table(kcut)
  out <- cbind(Count = count, Proportion = prop.table(count))
  structure(out, class = c("pareto_k_table", class(out)))
}

#' @rdname pareto-k-diagnostic
#' @export
#' @param digits The number of digits to print.
#' @param ... Ignored.
print.pareto_k_table <- function(x, digits = 1, ...) {
  count <- x[, "Count"]
  prop <- x[, "Proportion"]

  if (sum(count[2:4]) == 0) {
    cat("\nAll Pareto k estimates are good (k < 0.5)\n")
  } else {
    tab <- cbind(
      " " = rep("", 4),
      " " = c("(good)", "(ok)", "(bad)", "(very bad)"),
      "Count" = .fr(count, 0),
      " Pct" = paste0(.fr(100 * prop, digits), "%")
    )
    tab2 <- rbind(tab)
    cat("\nPareto k diagnostic values:\n")
    rownames(tab2) <- format(rownames(tab2), justify = "right")
    print(tab2, quote = FALSE)

    if (sum(count[3:4]) == 0)
      cat("\nAll Pareto k estimates are ok (k < 0.7)\n")

    cat(.k_help())
    invisible(x)
  }
}


# internal ----------------------------------------------------------------
get_pareto_k <- function(x) {
  if (is.null(x[["pareto_k"]]))
    stop("No Pareto k estimates found.", call. = FALSE)
  x[["pareto_k"]]
}
