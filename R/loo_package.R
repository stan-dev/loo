#' Efficient implementation of leave-one-out cross-validation and WAIC for
#' evaluating fitted Bayesian models
#'
#' @description Leave-one-out cross-validation and WAIC
#'
#' After fitting a Bayesian model we often want to measure its
#' predictive accuracy, for its own sake or for purposes of model comparison,
#' selection, or averaging (Geisser and Eddy, 1979, Hoeting et al., 1999,
#' Vehtari and Lampinen, 2002, Ando and Tsay, 2010, Vehtari and Ojanen, 2012).
#' Cross- validation and information criteria are two approaches for estimating
#' out-of-sample predictive accu- racy using within-sample fits (Akaike, 1973,
#' Stone, 1977). In this article we consider computations using the
#' log-likelihood evaluated at the usual posterior simulations of the
#' parameters. Computa- tion time for the predictive accuracy measures should be
#' negligible compared to the cost of fitting the model and obtaining posterior
#' draws in the first place. Exact cross-validation requires re-fitting the
#' model with different training sets. Approximate leave-one-out
#' cross-validation (LOO) can be computed easily using importance sampling
#' (Gelfand, Dey, and Chang, 1992, Gelfand, 1996) but the resulting estimate is
#' noisy, as the variance of the importance weights can be large or even
#' infinite (Peruggia, 1997, Epifani et al., 2008). Here we propose a novel
#' approach that provides a more accurate and reliable estimate using importance
#' weights that are smoothed using a Pareto distribution fit to the upper tail
#' of the distribution of importance weights. WAIC (the widely applicable or
#' Watanabe-Akaike information criterion; Watanabe, 2010) can be viewed as an
#' improvement on the deviance information criterion (DIC) for Bayesian models.
#' DIC has gained popularity in recent years in part through its implementation
#' in the graphical modeling package BUGS (Spiegelhalter, Best, et al., 2002;
#' Spiegelhalter, Thomas, et al., 1994, 2003), but it is known to have some
#' problems, arising in part from it not being fully Bayesian in that it is
#' based on a point estimate (van der Linde, 2005, Plummer, 2008). For example,
#' DIC can produce negative estimates of the effective number of parameters in a
#' model and it is not defined for singular models. WAIC is fully Bayesian and
#' closely approximates Bayesian cross-validation. Unlike DIC, WAIC is invariant
#' to parametrization and also works for singular models. WAIC is asymptotically
#' equal to LOO, and can thus be used as an approximation of LOO. In the finite
#' case, WAIC often gives similar estimates as LOO, but for influential
#' observations WAIC underestimates the effect of leaving out one observation.
#'
#' One advantage of AIC and DIC is their computational simplicity.
#' In the present paper, we quickly review LOO and WAIC and then present fast
#' and stable computations that can be per- formed directly on posterior
#' simulations, thus allowing these newer tools to enter routine statistical
#' practice. We compute LOO using very good importance sampling (VGIS), a new
#' procedure for regularizing importance weights (Vehtari and Gelman, 2015). As
#' a byproduct of our calculations, we also obtain approximate standard errors
#' for estimated predictive errors and for comparing of predictive errors
#' between two models.
#'
#' @docType package
#' @name loo
#'
NULL
