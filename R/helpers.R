#' More stable version of log(mean(exp(x)))
#'
#' @noRd
#' @param x A numeric vector.
#' @return A scalar equal to log(mean(exp(x))).
#'
logMeanExp <- function(x) {
  logS <- log(length(x))
  logSumExp(x) - logS
}

#' More stable version of log(colMeans(exp(x)))
#'
#' @noRd
#' @param x A matrix.
#' @return A vector where each element is LogSumExp of a column of x.
#'
colLogMeanExps <- function(x) {
  logS <- log(nrow(x))
  colLogSumExps(x) - logS
}

#' Compute point estimates and standard errors from pointwise vectors
#'
#' @noRd
#' @param x A matrix.
#' @return An ncol(x) by 2 matrix with columns 'Estimate' and 'SE'
#'   and rownames equal to colnames(x).
#'
table_of_estimates <- function(x) {
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

#' Check for NAs and non-finite values in log-lik (or log-ratios)
#' array/matrix/vector
#'
#' @noRd
#' @param x Array/matrix/vector of log-likelihood or log-ratio values.
#' @return x, invisibly, if no error is thrown.
#'
validate_ll <- function(x) {
  if (is.list(x)) {
    stop("List not allowed as input.")
  } else if (anyNA(x)) {
    stop("NAs not allowed in input.")
  } else if (!all(is.finite(x))) {
    stop("All input values must be finite.")
  }
  invisible(x)
}

#' Convert iter by chain by obs array to (iter * chain) by obs matrix
#'
#' @noRd
#' @param x Array to convert.
#' @return An (iter * chain) by obs matrix.
#'
llarray_to_matrix <- function(x) {
  stopifnot(is.array(x), length(dim(x)) == 3)
  xdim <- dim(x)
  dim(x) <- c(prod(xdim[1:2]), xdim[3])
  unname(x)
}

#' Convert (iter * chain) by obs matrix to iter by chain by obs array
#'
#' @noRd
#' @param x matrix to convert.
#' @param chain_id vector of chain ids.
#' @return iter by chain by obs array
#'
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


#' Validate that log-lik function exists and has correct arg names
#'
#' @noRd
#' @param x A function with arguments 'data_i' and 'draws'.
#' @return Either returns x or throws an error.
#'
validate_llfun <- function(x) {
  f <- match.fun(x)
  must_have <- c("data_i", "draws")
  arg_names <- names(formals(f))
  if (!all(must_have %in% arg_names)) {
    stop(
      "Log-likelihood function must have at least the arguments ",
      "'data_i' and 'draws'",
      call. = FALSE
    )
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


# Check how many cores to use and throw deprecation warning if loo.cores is used
loo_cores <- function(cores) {
  loo_cores_op <- getOption("loo.cores", NA)
  if (!is.na(loo_cores_op) && (loo_cores_op != cores)) {
    cores <- loo_cores_op
    warning("'loo.cores' is deprecated, please use 'mc.cores' or pass 'cores' explicitly.",
            call. = FALSE)
  }
  return(cores)
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
