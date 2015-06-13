#' Compare 
#' @export
#' @param loo1,loo2 lists returned from \code{\link{loo_and_waic}}. 
#' 

loo_and_waic_diff <- function(loo1, loo2){
  N1 <- nrow(loo1$pointwise)
  N2 <- nrow(loo2$pointwise)
  if (N1 != N2) {
    stop(paste("Models being compared should have the same number of data points.",
               "Found N1 =", N1, "and N2 =", N2))
  }
  sqrtN <- sqrt(N1)
  loo_diff <- loo2$pointwise[,"elpd_loo"] - loo1$pointwise[,"elpd_loo"]
  waic_diff <- loo2$pointwise[,"elpd_waic"] - loo1$pointwise[,"elpd_waic"]
  
  list(elpd_loo_diff = sum(loo_diff), lpd_loo_diff = sqrtN * sd(loo_diff),
       elpd_waic_diff = sum(waic_diff), lpd_waic_diff = sqrtN * sd(waic_diff))
}