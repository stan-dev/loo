#' @param is_method The importance sampling method to use.
#' The following methods are implemented:
#' * [`"psis"`][psis]: Pareto-Smoothed Importance Sampling (PSIS). Default method.
#' * [`"tis"`][tis]: Truncated Importance Sampling (TIS) with truncation at
#'   `sqrt(S)`, where `S` is the number of posterior draws.
#' * [`"sis"`][sis]: Standard Importance Sampling (SIS).
#'
