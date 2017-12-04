#' Helper functions for K-fold cross-validation
#'
#' These functions can be used to generate indexes for use with K-fold
#' cross-validation.
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
#' @details
#' \code{kfold_split_random} splits the data into \code{K} groups
#' of equal size (or roughly equal size).
#'
#' For a binary variable \code{x} that has many more \code{0}s than \code{1}s
#' (or vice-versa) \code{kfold_split_balanced} first splits the data by value of
#' \code{x}, does \code{kfold_split_random} within each of the two groups, and
#' then recombines the indexes returned from the two calls to
#' \code{kfold_split_random}. This helps ensure that the observations in the
#' less common category of \code{x} are more evenly represented across the
#' folds.
#'
#'
#'
#' @examples
#' kfold_split_random(K = 5, N = 20)
#'
#' y <- sample(c(0, 1), size = 200, replace = TRUE, prob = c(0.05, 0.95))
#' table(y)
#' ids <- kfold_split_balanced(K = 5, x = y)
#' table(ids[y == 0])
#'
#' grp <- gl(n = 50, k = 15, labels = state.name)
#' ids <- kfold_split_stratified(K = 10, grp)
#' cbind(ids, grp)
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
  S1 <- ceiling(Nlev / K) # Num levs in largest groups of levs
  N_S2 <- S1 * K - Nlev # Number of groups of levels of size S1 - 1
  N_S1 <- K - N_S2 # Number of groups of levels of size S1

  perm <- sample.int(Nlev) # permute group levels
  brks <- seq(from = S1 + 0.5, by = S1, length.out = N_S1)
  if (N_S2 > 0) {
    brks2 <- seq(from = brks[N_S1] + S1 - 1, by = S1 - 1, length.out = N_S2 - 1)
    brks <- c(brks, brks2)
  }
  grps <- findInterval(perm, vec = brks) + 1 # +1 so min is 1 not 0


  bins <- rep(NA, length(x))
  for (j in perm) {
    bins[x == j] <- grps[j]
  }
  return(bins)
}
