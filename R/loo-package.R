#' Efficient LOO-CV and WAIC for Bayesian models
#'
#' @docType package
#' @name loo-package
#'
#' @importFrom stats sd var quantile setNames weights rnorm qnorm
#' @importFrom matrixStats logSumExp colLogSumExps colSums2 colVars colMaxs
#'
#' @description
#' \if{html}{
#'   \figure{stanlogo.png}{options: width="50px" alt="mc-stan.org"}
#' }
#' *Stan Development Team*
#'
#' This package implements the methods described in Vehtari, Gelman, and
#' Gabry (2017), Vehtari, Simpson, Gelman, Yao, and Gabry (2019), and
#' Yao et al. (2018). To get started see the **loo** package
#' [vignettes](https://mc-stan.org/loo/articles/index.html), the
#' [loo()] function for efficient approximate leave-one-out
#' cross-validation (LOO-CV), the [psis()] function for the Pareto
#' smoothed importance sampling (PSIS) algorithm, or
#' [loo_model_weights()] for an implementation of Bayesian stacking of
#' predictive distributions from multiple models.
#'
#'
#' @details Leave-one-out cross-validation (LOO-CV) and the widely applicable
#'   information criterion (WAIC) are methods for estimating pointwise
#'   out-of-sample prediction accuracy from a fitted Bayesian model using the
#'   log-likelihood evaluated at the posterior simulations of the parameter
#'   values. LOO-CV and WAIC have various advantages over simpler estimates of
#'   predictive error such as AIC and DIC but are less used in practice because
#'   they involve additional computational steps. This package implements the
#'   fast and stable computations for approximate LOO-CV laid out in Vehtari,
#'   Gelman, and Gabry (2017a). From existing posterior simulation draws, we
#'   compute LOO-CV using Pareto smoothed importance sampling (PSIS; Vehtari,
#'   Simpson, Gelman, Yao, and Gabry, 2019), a new procedure for stabilizing
#'   and diagnosing importance weights. As a byproduct of our calculations,
#'   we also obtain approximate standard errors for estimated predictive
#'   errors and for comparing of predictive errors between two models.
#'
#'   We recommend PSIS-LOO-CV instead of WAIC, because PSIS provides useful
#'   diagnostics and effective sample size and Monte Carlo standard error
#'   estimates.
#'
#'
#' @template loo-and-psis-references
#' @template stacking-references
#' @template loo-large-data-references
#'
#' @references
#' Epifani, I., MacEachern, S. N., and Peruggia, M. (2008). Case-deletion
#' importance sampling estimators: Central limit theorems and related results.
#' *Electronic Journal of Statistics* **2**, 774-806.
#'
#' Gelfand, A. E. (1996). Model determination using sampling-based methods. In
#' *Markov Chain Monte Carlo in Practice*, ed. W. R. Gilks, S. Richardson,
#' D. J. Spiegelhalter, 145-162. London: Chapman and Hall.
#'
#' Gelfand, A. E., Dey, D. K., and Chang, H. (1992). Model determination using
#' predictive distributions with implementation via sampling-based methods. In
#' *Bayesian Statistics 4*, ed. J. M. Bernardo, J. O. Berger, A. P. Dawid,
#' and A. F. M. Smith, 147-167. Oxford University Press.
#'
#' Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. *Statistics and Computing*
#' **24**, 997-1016.
#'
#' Ionides, E. L. (2008). Truncated importance sampling. *Journal of
#' Computational and Graphical Statistics* **17**, 295-311.
#'
#' Koopman, S. J., Shephard, N., and Creal, D. (2009). Testing the assumptions
#' behind importance sampling. *Journal of Econometrics* **149**, 2-11.
#'
#' Peruggia, M. (1997). On the variability of case-deletion importance sampling
#' weights in the Bayesian linear model. *Journal of the American
#' Statistical Association* **92**, 199-207.
#'
#' Stan Development Team (2017). The Stan C++ Library, Version 2.17.0.
#' <https://mc-stan.org>.
#'
#' Stan Development Team (2018). RStan: the R interface to Stan, Version 2.17.3.
#' <https://mc-stan.org>.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely application information criterion in singular learning theory.
#' *Journal of Machine Learning Research* **11**, 3571-3594.
#'
#' Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
#' for the generalized Pareto distribution. *Technometrics* **51**,
#' 316-325.
#'
NULL
