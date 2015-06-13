#' Convenience function for extracting log-likelihood from 
#' a \code{stanfit} object. 
#' 
#' @export
#' @param stanfit a \code{stanfit} (\pkg{rstan}) object. 
#' @param parameter_name a character string naming the parameter (generated 
#' quantity) in the Stan model corresponding to the log-likelihood. 
#' @seealso \code{\link[rstan]{stanfit-class}}
log_lik <- function(stanfit, parameter_name = "log_lik") {
  rstan_ok <- requireNamespace("rstan", quietly = TRUE)
  if (!rstan_ok) {
    message("Please install the rstan package")
    return(invisible(NULL))
  }
  rstan::extract(stanfit, parameter_name)[[parameter_name]]
}

xx <- requireNamespace("rstan", quietly = TRUE)