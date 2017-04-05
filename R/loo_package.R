#' Efficient LOO and WAIC for Bayesian models
#'
#' @docType package
#' @name loo-package
#'
#' @importFrom stats sd var quantile
#'
#' @description
#' \if{html}{
#'   \figure{stanlogo.png}{options: width="50px" alt="mc-stan.org"}
#'   \emph{Stan Development Team}
#' }
#'
#' This package implements the methods described in Vehtari, Gelman,
#' and Gabry (2016a). The package documentation is largely based on excerpts
#' from the paper.
#'
#' @section Summary: Leave-one-out cross-validation (LOO) and the widely
#'   applicable information criterion (WAIC) are methods for estimating
#'   pointwise out-of-sample prediction accuracy from a fitted Bayesian model
#'   using the log-likelihood evaluated at the posterior simulations of the
#'   parameter values. LOO and WAIC have various advantages over simpler
#'   estimates of predictive error such as AIC and DIC but are less used in
#'   practice because they involve additional computational steps. This package
#'   implements the fast and stable computations for LOO and WAIC laid out in
#'   Vehtari, Gelman, and Gabry (2016a). From existing posterior simulation
#'   draws, we compute LOO using Pareto Smoothed Importance Sampling (PSIS;
#'   Vehtari, Gelman, and Gabry, 2016b), a new procedure for regularizing
#'   importance weights. As a byproduct of our calculations, we also obtain
#'   approximate standard errors for estimated predictive errors and for
#'   comparing of predictive errors between two models.
#'
#' @section Details: After fitting a Bayesian model we often want to measure its
#'   predictive accuracy, for its own sake or for purposes of model comparison,
#'   selection, or averaging. Cross-validation and information criteria are two
#'   approaches for estimating out-of-sample predictive accuracy using
#'   within-sample fits.
#'
#'   Exact cross-validation requires re-fitting the model with different
#'   training sets. Approximate leave-one-out cross-validation (LOO) can be
#'   computed easily using importance sampling (Gelfand, Dey, and Chang, 1992,
#'   Gelfand, 1996) but the resulting estimate is noisy, as the variance of the
#'   importance weights can be large or even infinite (Peruggia, 1997, Epifani
#'   et al., 2008). Here we propose a novel approach that provides a more
#'   accurate and reliable estimate using importance weights that are smoothed
#'   by fitting a generalized Pareto distribution to the upper tail of the
#'   distribution of the importance weights.
#'
#'   WAIC (the widely applicable or Watanabe-Akaike information criterion;
#'   Watanabe, 2010) can be viewed as an improvement on the deviance information
#'   criterion (DIC) for Bayesian models. DIC has gained popularity in recent
#'   years in part through its implementation in the graphical modeling package
#'   BUGS (Spiegelhalter, Best, et al., 2002; Spiegelhalter, Thomas, et al.,
#'   1994, 2003), but it is known to have some problems, arising in part from it
#'   not being fully Bayesian in that it is based on a point estimate (van der
#'   Linde, 2005, Plummer, 2008). For example, DIC can produce negative
#'   estimates of the effective number of parameters in a model and it is not
#'   defined for singular models. WAIC is fully Bayesian and closely
#'   approximates Bayesian cross-validation. Unlike DIC, WAIC is invariant to
#'   parametrization and also works for singular models. WAIC is asymptotically
#'   equal to LOO, and can thus be used as an approximation to LOO. In the
#'   finite case, WAIC and LOO often give very similar estimates, but for
#'   influential observations WAIC underestimates the effect of leaving out one
#'   observation.
#'
#'   One advantage of AIC and DIC has been their computational simplicity. In
#'   this package we present fast and stable computations for LOO and WAIC that
#'   can be performed directly on posterior simulations, thus allowing these
#'   newer tools to enter routine statistical practice. As a byproduct of our
#'   calculations, we also obtain approximate standard errors for estimated
#'   predictive errors and for the comparison of predictive errors between two
#'   models.
#'
#' @section PSIS-LOO: The distribution of the importance weights used in LOO may
#'   have a long right tail. We use the empirical Bayes estimate of Zhang and
#'   Stephens (2009) to fit a generalized Pareto distribution (gPd) to the tail
#'   (20\% largest importance ratios). By examining the shape parameter \eqn{k}
#'   of the fitted gPd, we are able to obtain sample based estimates of the
#'   existance of the moments (Koopman et al, 2009). This extends the diagnostic
#'   approach of Peruggia (1997) and Epifani et al. (2008) to be used routinely
#'   with importance-sampling LOO for any model with a factorizing likelihood.
#'
#'   Epifani et al. (2008) show that when estimating the leave-one-out
#'   predictive density, the central limit theorem holds if the variance of the
#'   weight distribution is finite. These results can be extended using the
#'   generalized central limit theorem for stable distributions. Thus, even if
#'   the variance of the importance weight distribution is infinite, if the mean
#'   exists the estimate's accuracy improves when additional draws are obtained.
#'   When the tail of the weight distribution is long, a direct use of
#'   importance sampling is sensitive to one or few largest values. By fitting a
#'   gPd to the upper tail of the importance weights we smooth these values. The
#'   procedure (implemented in the \code{\link{psislw}} function) goes as
#'   follows:
#'
#' \enumerate{
#'   \item Fit the gPd to the 20\% largest importance ratios \eqn{r_s}. The
#'   computation is done separately for each held-out data point \eqn{i}. In
#'   simulation experiments with thousands and tens of thousands of draws, we
#'   have found that the fit is not sensitive to the specific cutoff value (for
#'   a consistent estimation the proportion of the samples above the cutoff
#'   should get smaller when the number of draws increases).
#'
#'   \item Stabilize the importance ratios by replacing the \eqn{M} largest
#'   ratios by the expected values of the order statistics of the fitted
#'   gPd \deqn{G((z - 0.5)/M), z = 1,...,M,} where
#'   \eqn{M} is the number of simulation draws used to fit the Pareto (in this
#'   case, \eqn{M = 0.2*S}) and \eqn{G} is the inverse-CDF of the gPd.
#'
#'   \item To guarantee finite variance of the estimate, truncate the smoothed
#'   ratios with \deqn{S^{3/4}\bar{w},} where \eqn{\bar{w}} is the average of
#'   the smoothed weights.
#'}
#'
#' The above steps must be performed for each data point \eqn{i}. The result is
#' a vector of weights \eqn{w_{i}^{s}, s = 1,...,S}, for each \eqn{i}, which in
#' general should be better behaved than the raw importance ratios
#' \eqn{r_{i}^{s}} from which they were constructed.
#'
#' The results can be then combined to compute the desired LOO estimates.
#'
#' @section Pareto k diagnostic: The reliability of the PSIS-based estimates can
#'   be assessed using the estimates for the shape parameter \eqn{k} of the
#'   generalized Pareto distribution.
#'
#' \itemize{
#'   \item If \eqn{k < 1/2} the variance of the raw importance ratios is finite,
#'   the central limit theorem holds, and the estimate converges quickly.
#'
#'   \item If \eqn{k} is between 1/2 and 1 the variance of the raw importance
#'   ratios is infinite but the mean exists, the generalized central limit
#'   theorem for stable distributions holds, and the convergence of the estimate
#'   is slower. The variance of the PSIS estimate is finite but may be large.
#'
#'   \item If \eqn{k > 1} the variance and the mean of the raw ratios
#'   distribution do not exist. The variance of the PSIS estimate is finite but
#'   may be large.
#' }
#'
#' If the estimated tail shape parameter \eqn{k} exceeds \eqn{0.5}, the user
#' should be warned, although in practice we have observed good performance for
#' values of \eqn{k} up to 0.7. Even if the PSIS estimate has a finite variance,
#' the user should consider sampling directly from \eqn{p(\theta^s | y_{-i})}
#' for the problematic \eqn{i}, use \eqn{k}-fold cross-validation, or use a more
#' robust model.
#'
#' Importance sampling is likely to work less well if the marginal posterior
#' \eqn{p(\theta^s | y)} and LOO posterior \eqn{p(\theta^s | y_{-i})} are much
#' different, which is more likely to happen with a non-robust model and highly
#' influential observations. A robust model may reduce the sensitivity to highly
#' influential observations.
#'
#' @template loo-and-psis-references
#' @references
#' Epifani, I., MacEachern, S. N., and Peruggia, M. (2008). Case-deletion
#' importance sampling estimators: Central limit theorems and related results.
#' \emph{Electronic Journal of Statistics} \strong{2}, 774-806.
#'
#' Gelfand, A. E. (1996). Model determination using sampling-based methods. In
#' \emph{Markov Chain Monte Carlo in Practice}, ed. W. R. Gilks, S. Richardson,
#' D. J. Spiegelhalter, 145-162. London: Chapman and Hall.
#'
#' Gelfand, A. E., Dey, D. K., and Chang, H. (1992). Model determination using
#' predictive distributions with implementation via sampling-based methods. In
#' \emph{Bayesian Statistics 4}, ed. J. M. Bernardo, J. O. Berger, A. P. Dawid,
#' and A. F. M. Smith, 147-167. Oxford University Press.
#'
#' Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}
#' \strong{24}, 997-1016.
#'
#' Ionides, E. L. (2008). Truncated importance sampling. \emph{Journal of
#' Computational and Graphical Statistics} \strong{17}, 295-311.
#'
#' Koopman, S. J., Shephard, N., and Creal, D. (2009). Testing the assumptions
#' behind importance sampling. \emph{Journal of Econometrics} \strong{149}, 2-11.
#'
#' Peruggia, M. (1997). On the variability of case-deletion importance sampling
#' weights in the Bayesian linear model. \emph{Journal of the American
#' Statistical Association} \strong{92}, 199-207.
#'
#' Stan Development Team (2016). The Stan C++ Library, Version 2.10.0.
#' \url{http://mc-stan.org/documentation/}.
#'
#' Stan Development Team (2016). RStan: the R interface to Stan, Version 2.10.1
#' \url{http://mc-stan.org/interfaces/rstan.html}.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely application information criterion in singular learning theory.
#' \emph{Journal of Machine Learning Research} \strong{11}, 3571-3594.
#'
#' Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
#' for the generalized Pareto distribution. \emph{Technometrics} \strong{51},
#' 316-325.
#'
NULL
