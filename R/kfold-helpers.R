#' Helper functions for K-fold cross-validation
#'
#' @name kfold-helpers
#' @param K The number of folds to use.
#' @param N The number of observations in the data.
#' @param x A discrete variable of length \code{N}. Will be coerced to
#'   \code{\link{factor}}. For \code{kfold_split_balanced} \code{x} should be a
#'   binary variable. For \code{kfold_split_stratified} \code{x} should be a
#'   grouping variable with at least \code{K} levels.
#' @return An integer vector of length \code{N} where each element is an index
#'   in \code{1:K}.
#'
NULL

#' @rdname kfold-helpers
#' @export
kfold_split_random <- function(K = 10, N = NULL) {
  perm <- sample.int(N)
  idx <- ceiling(seq(from = 1, to = N, length.out = K + 1))
  bin <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
  return(bin)
}

#' @rdname kfold-helpers
#' @export
kfold_split_balanced <- function(K = 10, x = NULL) {
  stopifnot(length(unique(x)) == 2, K <= length(x))
  x <- as.integer(as.factor(x)) - 1
  x0_id <- which(x == 0)
  x1_id <- which(x == 1)
  bins <- rep(NA, length(x))
  bins[x0_id] <- kfold_split_random(K = K, N = length(x0_id))
  bins[x1_id] <- kfold_split_random(K = K, N = length(x1_id))
  return(bins)
}

#' @rdname kfold-helpers
#' @export
kfold_split_stratified <- function(K = 10, x = NULL) {
  Nlev <- length(unique(x))
  if (Nlev < K) {
    stop("K must be at least the number of groups in x.")
  }
  x <- as.integer(as.factor(x))
  if (Nlev == K) {
    return(x)
  }

  # Otherwise we have Nlev > K
  Size1 <- ceiling(Nlev / K) # Num levs in largest groups of levs
  N_Size2 <- Size1 * K - Nlev # Number of groups of levels of size Size1 - 1
  N_Size1 <- K - N_Size2 # Number of groups of levels of size S1
  brks <- rep(NA, K-1)
  brks[1] <- Size1 + 0.5
  for (j in 2:N_Size1) {
    brks[j] <- brks[j-1] + Size1
  }
  for (j in (N_Size1 + 1):(K-1)) {
    brks[j] <- brks[j-1] + (Size1 - 1)
  }

  perm <- sample.int(Nlev) # permute group levels
  grps <- findInterval(perm, brks) + 1 # +1 so min is 1 not 0
  bins <- rep(NA, length(x))
  for (j in seq_along(perm)) {
    bins[perm == j] <- grps[j]
  }
  return(bins)
}
