#' @param is_method
#'   Importance sampling method to use. The following methods are implemented:
#' \describe{
#'  \item{\code{psis}}{
#'    Pareto-Smoothed Importance Sampling (PSIS). Default method.
#'  }
#'  \item{\code{tis}}{
#'    Truncated Importance Sampling (TIS) with truncation at \code{sqrt(S)},
#'     where \code{S} is the number of posterior draws.
#'  }
#'  \item{\code{sis}}{
#'    Standard Importance Sampling (SIS).
#'  }
#' }
#'

