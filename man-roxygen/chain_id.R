#' @param chain_id A vector of length \code{NROW(x)} containing MCMC chain
#'   indexes for each each row of \code{x} (if a matrix) or each value in
#'   \code{x} (if a vector). No \code{chain_id} is needed if \code{x} is a 3-D
#'   array. If there are \code{C} chains then valid chain indexes are values
#'   in \code{1:C}.
