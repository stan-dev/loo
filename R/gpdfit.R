#' Estimate parameters of the Generalized Pareto distribution
#'
#' Given a sample \eqn{x}, Estimate the parameters \eqn{k} and \eqn{\sigma} of
#' the generalized Pareto distribution (GPD), assuming the location parameter is
#' 0. By default the fit uses a prior for \eqn{k}, which will stabilize
#' estimates for very small sample sizes (and low effective sample sizes in the
#' case of MCMC samples). The weakly informative prior is a Gaussian prior
#' centered at 0.5.
#'
#' @export
#' @param x A numeric vector. The sample from which to estimate the parameters.
#' @param wip Logical indicating whether to adjust \eqn{k} based on a weakly
#'   informative Gaussian prior centered on 0.5. Defaults to `TRUE`.
#' @param min_grid_pts The minimum number of grid points used in the fitting
#'   algorithm. The actual number used is `min_grid_pts + floor(sqrt(length(x)))`.
#' @param sort_x If `TRUE` (the default), the first step in the fitting
#'   algorithm is to sort the elements of `x`. If `x` is already
#'   sorted in ascending order then `sort_x` can be set to `FALSE` to
#'   skip the initial sorting step.
#' @return A named list with components `k` and `sigma`.
#'
#' @details Here the parameter \eqn{k} is the negative of \eqn{k} in Zhang &
#'   Stephens (2009).
#'
#' @seealso [psis()], [pareto-k-diagnostic]
#'
#' @references
#' Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
#' for the generalized Pareto distribution. *Technometrics* **51**, 316-325.
#'
gpdfit <- function(x, wip = TRUE, min_grid_pts = 30, sort_x = TRUE) {
  posterior::gpdfit(
    x = x,
    wip = wip,
    min_grid_pts = min_grid_pts,
    sort_x = sort_x
  )
}
