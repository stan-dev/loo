#' Objects to use in examples and tests
#'
#' Example pointwise log-likelihood objects to use in demonstrations and tests.
#' See the \strong{Value} and \strong{Examples} sections below.
#'
#' @export
#'
#' @return
#' \code{example_loglik_array} returns a 500 (draws) x 2 (chains) x 32
#' (observations) pointwise log-likelihood array.
#'
#' \code{example_loglik_matrix} returns the same pointwise log-likelihood values
#' as \code{example_loglik_array} but reshaped into a 1000 (draws*chains) x 32
#' (observations) matrix.
#'
#' @examples
#' LLarr <- example_loglik_array()
#' (dim_arr <- dim(LLarr))
#' LLmat <- example_loglik_matrix()
#' (dim_mat <- dim(LLmat))
#'
#' all.equal(dim_mat[1], dim_arr[1] * dim_arr[2])
#' all.equal(dim_mat[2], dim_arr[3])
#'
#' all.equal(LLarr[, 1, ], LLmat[1:500, ])
#' all.equal(LLarr[, 2, ], LLmat[501:1000, ])
#'
example_loglik_array <- function() {
  # .example_loglik_array exists in R/sysdata.R
  return(.example_loglik_array)
}

#' @rdname example_loglik_array
#' @export
example_loglik_matrix <- function() {
  ll <- example_loglik_array()
  return(llarray_to_matrix(ll))
}

