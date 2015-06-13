#' Very good importance sampling
#' 
#' @export
#' @param log_lik, wcp, wtrunc arguments.
#' @return a list.
#' @details The distribution of the importance weights used in LOO may have a
#'   long right tail. We use the empirical Bayes estimate of Zhang and Stephens
#'   (2009) to fit a generalized Pareto distribution to the tail (20% largest
#'   importance ratios). By examining the shape parameter k of the fitted Pareto
#'   distribution, we are able to obtain sample based estimate of the existence
#'   of the moments (Koopman et al, 2009). This extends the diagnostic approach
#'   of Peruggia (1997) and Epifani et al. (2008) to be used routinely with
#'   IS-LOO for any model with factorising likelihood. Epifani et al. (2008)
#'   show that when estimating the leave-one-out predictive density, the central
#'   limit theorem holds if the variance of the weight distribution is finite.
#'   These results can be extended by using the generalized central limit
#'   theorem for stable distributions. Thus, even if the variance of the
#'   importance weight distribution is infinite, if the mean exists the
#'   estimate’s accuracy improves when additional draws are obtained. When the
#'   tail of the weight distribution is long, a direct use of importance
#'   sampling is sensitive to one or few largest values. By fitting a
#'   generalized Pareto distribution to the upper tail of the importance
#'   weights, we smooth these values. The procedure goes as follows:
#'   
#'   \enumerate{
#'   \item Fit the generalized Pareto distribution to the 20% largest importance
#'   ratios \eqn{r_s} as computed in (6). (The computation is done separately for each
#'   held-out data point \eqn{i}.) In simulation experiments with a thousands to tens
#'   of thousands of simulation draws, we have found the fit is not sensitive to
#'   the specific cutoff value (for a consistent estimation the proportion of
#'   the samples above the cutoff should get smaller when the number of draws
#'   increases).
#'   
#'   \item Stabilize the importance ratios by replacing the \eqn{M} largest ratios 
#'   by the expected values of the order statistics of the fitted generalized
#'   Pareto distribution \deqn{G((z - 0.5)/M), z = 1,...,M,}
#'   where \eqn{M} is the number of simulation draws used to fit the Pareto (in this
#'   case, \eqn{M = 0.2*S}) and \eqn{G} is the inverse-CDF of the generalized
#'   Pareto distribution.
#'   
#'   \item To guarantee finite variance of the estimate, truncate the smoothed
#'   ratios with \deqn{S^{3/4}\bar{w},} where \eqn{\bar{w}} is the average of 
#'   the smoothed weights.
#'   }
#'   
#'   The above steps must be performed for each data point \eqn{i}, thus 
#'   resulting in a vector of weights \eqn{w_{i}^{s}, s = 1,...,S}, for each 
#'   \eqn{i}, which in general should be better behaved than the raw importance 
#'   ratios \eqn{r_{i}^{s}} from which they were constructed.
#'   
#'   The results can then be combined to compute desired LOO estimates.
#'

vgisloo <- function(log_lik, wcp=20, wtrunc=3/4) {
  lw <- -log_lik
  temp <- vgislw(lw, wcp, wtrunc)
  vglw <- temp$lw
  vgk <- temp$k
  loos <- sumlogs(log_lik + vglw)
  loo <- sum(loos)
  list(loo=loo, loos=loos, ks=vgk)
}