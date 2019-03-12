#' Helper functions for K-fold cross-validation
#'
#' @description These functions can be used to generate indexes for use with
#'   K-fold cross-validation. See the \strong{Details} section for explanations.
#'
#' @name kfold-helpers
#' @param K The number of folds to use.
#' @param N The number of observations in the data.
#' @param x A discrete variable of length \code{N} with at least \code{K} levels
#'   (unique values). Will be coerced to \code{\link{factor}}.
#' 
#' @return An integer vector of length \code{N} where each element is an index in \code{1:K}.
#'
#' @details
#' \code{kfold_split_random} splits the data into \code{K} groups
#' of equal size (or roughly equal size).
#'
#' For a categorical variable \code{x} \code{kfold_split_stratified}
#' splits the observations into \code{K} groups ensuring that relative
#' category frequencies are approximately preserved.
#'
#' For a grouping variable \code{x}, \code{kfold_split_grouped} places
#' all observations in \code{x} from the same group/level together in
#' the same fold. The selection of which groups/levels go into which
#' fold (relevant when when there are more folds than groups) is
#' randomized.
#'
#' @examples
#' ids <- kfold_split_random(K = 5, N = 20)
#' print(ids)
#' table(ids)
#'
#'
#' x <- sample(c(0, 1), size = 200, replace = TRUE, prob = c(0.05, 0.95))
#' table(x)
#' ids <- kfold_split_stratified(K = 5, x = x)
#' print(ids)
#' table(ids, x)
#'
#' grp <- gl(n = 50, k = 15, labels = state.name)
#' length(grp)
#' head(table(grp))
#'
#' ids_10 <- kfold_split_grouped(K = 10, x = grp)
#' (tab_10 <- table(grp, ids_10))
#' colSums(tab_10)
#'
#' ids_9 <- kfold_split_grouped(K = 9, x = grp)
#' (tab_9 <- table(grp, ids_9))
#' colSums(tab_9)
#'
NULL

#' @rdname kfold-helpers
#' @export
kfold_split_random <- function(K = 10, N = NULL) {
  stopifnot(
    !is.null(N),
    K == as.integer(K),
    N == as.integer(N),
    length(K) == 1,
    length(N) == 1,
    K > 1,
    K <= N
  )
  perm <- sample.int(N)
  idx <- ceiling(seq(from = 1, to = N, length.out = K + 1))
  bins <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
  return(bins)
}

#' @rdname kfold-helpers
#' @export
kfold_split_stratified <- function(K = 10, x = NULL) {
  stopifnot(
    !is.null(x),
    K == as.integer(K),
    length(K) == 1,
    K > 1,
    K <= length(x)
  )
  x <- as.integer(as.factor(x))
  Nlev <- length(unique(x))
  N <- length(x)
  xids <- numeric()
  for (l in 1:Nlev) {
      xids <- c(xids, sample(which(x==l)))
  }
  bins <- rep(NA, N)
  bins[xids] <- rep(1:K, ceiling(N/K))[1:N]
  return(bins)
}

#' @rdname kfold-helpers
#' @export
kfold_split_grouped <- function(K = 10, x = NULL) {
  stopifnot(
    !is.null(x),
    K == as.integer(K),
    length(K) == 1,
    K > 1,
    K <= length(x)
  )

  Nlev <- length(unique(x))
  if (Nlev < K) {
    stop("'K' must not be bigger than the number of levels/groups in 'x'.")
  }
  x <- as.integer(as.factor(x))
  if (Nlev == K) {
    return(x)
  }

  # Otherwise we have Nlev > K
  S1 <- ceiling(Nlev / K)  # number of levels in largest groups of levels
  N_S2 <- S1 * K - Nlev    # number of groups of levels of size S1 - 1
  N_S1 <- K - N_S2         # number of groups of levels of size S1

  perm <- sample.int(Nlev) # permute group levels
  brks <- seq(from = S1 + 0.5, by = S1, length.out = N_S1)
  if (N_S2 > 0) {
    brks2 <- seq(from = brks[N_S1] + S1 - 1, by = S1 - 1, length.out = N_S2 - 1)
    brks <- c(brks, brks2)
  }
  grps <- findInterval(perm, vec = brks) + 1  # +1 so min is 1 not 0

  bins <- rep(NA, length(x))
  for (j in perm) {
    bins[x == j] <- grps[j]
  }
  return(bins)
}
