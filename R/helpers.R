#' More stable version of `log(mean(exp(x)))`
#'
#' @noRd
#' @param x A numeric vector.
#' @return A scalar equal to `log(mean(exp(x)))`.
#'
logMeanExp <- function(x) {
  logS <- log(length(x))
  matrixStats::logSumExp(x) - logS
}

#' More stable version of `log(colMeans(exp(x)))`
#'
#' @noRd
#' @param x A matrix.
#' @return A vector where each element is `logMeanExp()` of a column of `x`.
#'
colLogMeanExps <- function(x) {
  logS <- log(nrow(x))
  matrixStats::colLogSumExps(x) - logS
}

#' More stable version of `exp(x) - exp(y)`
#'
#' @noRd
#' @param x A numeric vector.
#' @param y A numeric scalar or vector recycled to the length of `x`.
#'   Must satisfy `x >= y` elementwise.
#' @return A numeric vector equal to `exp(x) - exp(y)`.
#'
exp_x_minus_exp_y <- function(x, y) {
  out <- -exp(x) * expm1(y - x)
  # which() drops the NA comparisons that NA or NaN inputs would produce
  out[which(x == y)] <- 0
  out
}

#' More stable version of `x^2 - y^2`
#'
#' @noRd
#' @param x,y Numeric vectors of the same length.
#' @return A numeric vector equal to `x^2 - y^2`.
#'
difference_of_squares <- function(x, y) {
  (x - y) * (x + y)
}

#' More stable version of `(exp(a) - exp(b)) / exp(c)`
#'
#' @noRd
#' @param a,b,c Numeric vectors of the same length.
#' @return A numeric vector equal to `(exp(a) - exp(b)) / exp(c)`. Elements
#'   with `a == b` are returned as an exact zero regardless of `c`; elsewhere
#'   `NA` and `NaN` inputs propagate.
#'
exp_diff_over_exp <- function(a, b, c) {
  # `a >= b` is NA if `a` or `b` is NA or NaN, and R silently ignores NA
  # indices in `[<-`. Seed the result from the inputs and index with which()
  # so that missing values propagate instead of leaving a zero behind.
  out <- a + b + c
  larger <- which(a >= b)
  smaller <- which(a < b)
  out[larger] <- exp(a[larger] - c[larger]) * -expm1(b[larger] - a[larger])
  out[smaller] <- exp(b[smaller] - c[smaller]) * expm1(a[smaller] - b[smaller])
  out[which(a == b)] <- 0
  out
}

#' Compute point estimates and standard errors from pointwise vectors
#'
#' @noRd
#' @param x A matrix.
#' @return An `ncol(x)` by 2 matrix with columns `"Estimate"` and `"SE"`
#'   and rownames equal to `colnames(x)`.
#'
table_of_estimates <- function(x) {
  out <- cbind(
    Estimate = matrixStats::colSums2(x),
    SE = sqrt(nrow(x) * matrixStats::colVars(x))
  )
  rownames(out) <- colnames(x)
  return(out)
}


# validating and reshaping arrays/matrices  -------------------------------

#' Check for `NA` and non-finite values in log-lik (or log-ratios)
#' array/matrix/vector
#'
#' @noRd
#' @param x Array/matrix/vector of log-likelihood or log-ratio values.
#' @return `x`, invisibly, if no error is thrown.
#'
validate_ll <- function(x) {
  if (is.list(x)) {
    stop("List not allowed as input.")
  } else if (anyNA(x)) {
    stop("NAs not allowed in input.")
  } else if (any(x == Inf)) {
    stop("All input values must be finite or -Inf.")
  }
  invisible(x)
}

#' Check that a log-likelihood array/matrix/vector is finite
#'
#' `loo()` requires finite log-likelihood values, for two unrelated reasons.
#'
#' A `-Inf` log likelihood is meaningful on its own — the observation has zero
#' likelihood under that draw — but the leave-one-out importance ratio is
#' `1 / p(y_i | theta)`, which is then infinite, so the PSIS estimate does not
#' exist. Because `loo()` negates the log likelihood before importance
#' sampling, this used to surface as [validate_ll()]'s `+Inf` message, with the
#' polarity reversed.
#'
#' A `+Inf` log likelihood is not meaningful, and it passes the log-ratio check
#' as a `-Inf` ratio. `ll + lw` is then `Inf + -Inf`, so a single such value
#' made every estimate for the model `NA`.
#'
#' This is deliberately stricter than [validate_ll()], which is also used for
#' log ratios, where `-Inf` is a valid zero importance weight.
#'
#' @noRd
#' @param x Array/matrix/vector of log-likelihood values.
#' @return `x`, invisibly, if no error is thrown.
#'
validate_log_lik <- function(x) {
  if (is.list(x)) {
    stop("List not allowed as input.")
  }
  # single pass covering NA, NaN and both infinities; the more specific checks
  # below only run when something is already known to be wrong
  if (!all(is.finite(x))) {
    if (anyNA(x)) {
      stop("NAs not allowed in input.")
    }
    stop("All log-likelihood values must be finite.")
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
#' @param x A function with arguments `data_i` and `draws`.
#' @return Either returns `x` or throws an error.
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
#' names of the objects in the list. The code `list(a = a, b = b)` becomes
#' `nlist(a,b)` and `list(a = a, b = 2)` becomes `nlist(a, b = 2)`, etc.
#'
#' @export
#' @keywords internal
#' @param ... Objects to include in the list.
#' @return A named list.
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

#' Map a function over `X`, optionally on a mizu pool or mirai daemons
#'
#' With `parallel = FALSE` this is `lapply()`. With `parallel = TRUE`, `X` is
#' split into blocks that run on the default mizu pool if one is set with
#' [mizu::mizu_set_default_pool()], and otherwise on the daemons set up with
#' [mirai::daemons()]. The results come back in the original order.
#'
#' `FUN` and `...` are sent to the workers, so `FUN` should not capture a large
#' environment; pass data as arguments instead. mizu workers and local mirai
#' daemons read the arguments from shared memory ([mori::share()]); remote
#' daemons receive copies. Warnings raised on a worker are signalled again on
#' the host, and the first error stops the map.
#'
#' @noRd
#' @param X A vector or list to map over.
#' @param FUN The function to apply to each element of `X`.
#' @param ... Further arguments passed to `FUN`.
#' @param parallel If `TRUE`, run on the mizu pool or mirai daemons.
#' @param batch If `TRUE`, send a few blocks of consecutive elements per worker,
#'   which suits many cheap elements. If `FALSE`, send each element as its own
#'   task so that idle workers take the next one, which suits a few expensive,
#'   uneven elements.
#' @return A list the length of `X`, as from `lapply()`.
loo_map <- function(X, FUN, ..., parallel, batch = TRUE) {
  if (!parallel || length(X) == 0L) {
    return(lapply(X, FUN, ...))
  }
  blocks <- if (mizu_pool_is_set()) {
    map_blocks_mizu(X, FUN, ..., batch = batch)
  } else {
    map_blocks_mirai(X, FUN, ..., batch = batch)
  }
  lapply(concat_field(blocks, "warnings"), warning)
  concat_field(blocks, "value")
}

map_blocks_mirai <- function(X, FUN, ..., batch) {
  mirai::require_daemons()
  n_blocks <- if (batch) mirai::info()[["connections"]] else length(X)
  share <- if (daemons_are_local()) mori::share else identity
  mirai::mirai_map(
    split_blocks(X, n_blocks),
    lapply_with_warnings,
    .args = c(list(.fun = FUN), lapply(list(...), share))
  )[mirai::.stop]
}

map_blocks_mizu <- function(X, FUN, ..., batch) {
  pool <- mizu::mizu_default_pool()
  n_workers <- sum(mizu::mizu_pool_info(pool)[["workers"]] == "live")
  n_blocks <- if (batch) n_workers * mizu_blocks_per_worker else length(X)
  groups <- split_blocks(X, n_blocks)
  handle <- do.call(
    mizu::mizu_map,
    c(
      list(pool, groups, lapply_with_warnings, .fun = FUN, .n_chunks = length(groups)),
      lapply(list(...), mori::share)
    )
  )
  mizu::mizu_collect(handle)
}

mizu_blocks_per_worker <- 8L

mizu_pool_is_set <- function() {
  isNamespaceLoaded("mizu") && !is.null(mizu::mizu_default_pool())
}

split_blocks <- function(X, n_blocks) {
  n_blocks <- max(1L, min(n_blocks, length(X)))
  unname(split(X, ceiling(seq_along(X) * n_blocks / length(X))))
}

#' Split a matrix or 3-d array along its last dimension
#'
#' mori shares a list of vectors without copying it, but gives every daemon a
#' private copy of a whole array as soon as the array is indexed. Daemons
#' therefore receive `slices(x, parallel)` and read it with `slice(x, i)`.
#' Sequential calls keep `x` as it is.
#'
#' @noRd
slices <- function(x, parallel) {
  if (!parallel) {
    return(x)
  }
  lapply(seq_len(dim(x)[length(dim(x))]), slice, x = x)
}

slice <- function(x, i) {
  if (is.list(x)) x[[i]] else if (length(dim(x)) == 2L) x[, i] else x[, , i]
}

daemons_are_local <- function() {
  isTRUE(grepl("^(ipc|abstract)://", mirai::nextget("url")))
}

concat_field <- function(blocks, field) {
  unlist(lapply(blocks, `[[`, field), recursive = FALSE)
}

lapply_with_warnings <- function(block, .fun, ...) {
  warnings <- list()
  value <- withCallingHandlers(
    lapply(block, .fun, ...),
    warning = function(w) {
      warnings <<- c(warnings, list(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}


# nocov start
# release reminders (for devtools)
release_questions <- function() {
  c(
    "Have you updated references?",
    "Have you updated inst/CITATION?",
    "Have you updated the vignettes?"
  )
}
# nocov end

is_constant <- function(x, tol = .Machine$double.eps) {
  abs(max(x) - min(x)) < tol
}
