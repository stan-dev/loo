#' @importFrom matrixStats logSumExp colLogSumExps colSums2 colVars

# more stable version of log(colMeans(exp(x)))
# @param x matrix
logColMeansExp <- function(x) {
  logS <- log(nrow(x))
  colLogSumExps(x) - logS
}

# more stable version of log(mean(exp(x)))
logMeanExp <- function(x) {
  logSumExp(x) - log(length(x))
}

#' Compute point estimates and standard errors from pointwise vectors
#'
#' @noRd
#' @param x A matrix.
#' @return An ncol(pointwise) by 2 matrix with columns 'Estimate' and 'SE'
#'   and rownames equal to colnames(pointwise).
#'
table_of_estimates <- function(x) {
  if ("mcse_elpd_loo" %in% colnames(x)) {
    x <- x[, !colnames(x) %in% "mcse_elpd_loo"]
  }
  out <- cbind(
    Estimate = colSums2(x),
    SE = sqrt(nrow(x) * colVars(x))
  )
  rownames(out) <- colnames(x)
  return(out)
}


# checking classes --------------------------------------------------------
is.psis <- function(x) {
  inherits(x, "psis") && is.list(x)
}
is.loo <- function(x) {
  inherits(x, "loo")
}
is.psis_loo <- function(x) {
  inherits(x, "psis_loo") && is.loo(x)
}
is.waic <- function(x) {
  inherits(x, "waic") && is.loo(x)
}


# validating and reshaping arrays/matrices  -------------------------------

# check for NAs and non-finite values in log-lik (or log-weights) array/matrix/vector
# @param x log-lik array/matrix/vector
# @return x if no error is thrown.
#
validate_ll <- function(x) {
  stopifnot(!is.list(x), !anyNA(x), all(is.finite(x)))
  invisible(x)
}

# Convert iter by chain by obs array to (iter * chain) by obs matrix
#
# @param x array to convert.
# @return (iter * chain) by obs matrix
#
llarray_to_matrix <- function(x) {
  stopifnot(is.array(x), length(dim(x)) == 3)
  xdim <- dim(x)
  dim(x) <- c(prod(xdim[1:2]), xdim[3])
  unname(x)
}

# Convert (iter * chain) by obs matrix to iter by chain by obs array
#
# @param x matrix to convert.
# @param chain_id vector of chain ids.
# @return iter by chain by obs array
#
llmatrix_to_array <- function(x, chain_id) {
  stopifnot(is.matrix(x), all(chain_id == as.integer(chain_id)))

  lldim <- dim(x)
  n_chain <- length(unique(chain_id))
  chain_id <- as.integer(chain_id)
  chain_counts <- as.numeric(table(chain_id))

  if (length(chain_id) != lldim[1]) {
    stop("Number of rows in matrix not equal to length(chain_id).",
         call. = FALSE)
  } else if (any(chain_counts != chain_counts[1])) {
    stop("Not all chains have same number of iterations.",
         call. = FALSE)
  } else if (max(chain_id) != n_chain) {
    stop("max(chain_id) not equal to the number of chains.",
         call. = FALSE)
  }

  n_iter <- lldim[1] / n_chain
  n_obs <- lldim[2]
  a <- array(data = NA, dim = c(n_iter, n_chain, n_obs))
  for (c in seq_len(n_chain)) {
    a[, c, ] <- x[chain_id == c, , drop = FALSE]
  }
  return(a)
}


# validate log-lik function
validate_llfun <- function(x) {
  f <- match.fun(x)
  arg_names <- names(formals(f))
  if (length(arg_names) != 2 ||
      !all(arg_names %in% c("data_i", "draws"))) {
    stop("Log-likelihood function should have two arguments: ",
         "'data_i' and 'draws'.", call. = FALSE)
  }
  return(f)
}



#' Named lists
#'
#' Create a named list using specified names or, if names are omitted, using the
#' names of the objects in the list. The code \code{list(a = a, b = b)} becomes
#' \code{nlist(a,b)} and \code{list(a = a, b = 2)} becomes \code{nlist(a, b =
#' 2)}, etc.
#'
#' @export
#' @keywords internal
#' @param ... Objects to include in the list.
#' @return A named list.
#'
#' @seealso \code{\link[base]{list}}
#' @author Jonah Gabry
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
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }

  return(out)
}

# nocov start
# release reminders (for devtools)
release_questions <- function() {
  c(
    "Have you updated all references to the LOO paper?",
    "Have you updated inst/CITATION?",
    "Have you updated R code in vignette to match the code in the paper?"
  )
}
# nocov end
