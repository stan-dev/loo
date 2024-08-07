% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loo-package.R
\docType{package}
\name{loo-package}
\alias{loo-package}
\title{Efficient LOO-CV and WAIC for Bayesian models}
\description{
\if{html}{
  \figure{stanlogo.png}{options: width="50" alt="mc-stan.org"}
}
\emph{Stan Development Team}

This package implements the methods described in Vehtari, Gelman, and
Gabry (2017), Vehtari, Simpson, Gelman, Yao, and Gabry (2024), and
Yao et al. (2018). To get started see the \strong{loo} package
\href{https://mc-stan.org/loo/articles/index.html}{vignettes}, the
\code{\link[=loo]{loo()}} function for efficient approximate leave-one-out
cross-validation (LOO-CV), the \code{\link[=psis]{psis()}} function for the Pareto
smoothed importance sampling (PSIS) algorithm, or
\code{\link[=loo_model_weights]{loo_model_weights()}} for an implementation of Bayesian stacking of
predictive distributions from multiple models.
}
\details{
Leave-one-out cross-validation (LOO-CV) and the widely applicable
information criterion (WAIC) are methods for estimating pointwise
out-of-sample prediction accuracy from a fitted Bayesian model using the
log-likelihood evaluated at the posterior simulations of the parameter
values. LOO-CV and WAIC have various advantages over simpler estimates of
predictive error such as AIC and DIC but are less used in practice because
they involve additional computational steps. This package implements the
fast and stable computations for approximate LOO-CV laid out in Vehtari,
Gelman, and Gabry (2017). From existing posterior simulation draws, we
compute LOO-CV using Pareto smoothed importance sampling (PSIS; Vehtari,
Simpson, Gelman, Yao, and Gabry, 2024), a new procedure for stabilizing
and diagnosing importance weights. As a byproduct of our calculations,
we also obtain approximate standard errors for estimated predictive
errors and for comparing of predictive errors between two models.

We recommend PSIS-LOO-CV instead of WAIC, because PSIS provides useful
diagnostics and effective sample size and Monte Carlo standard error
estimates.
}
\references{
Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC.
\emph{Statistics and Computing}. 27(5), 1413--1432. doi:10.1007/s11222-016-9696-4
(\href{https://link.springer.com/article/10.1007/s11222-016-9696-4}{journal version},
\href{https://arxiv.org/abs/1507.04544}{preprint arXiv:1507.04544}).

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
Pareto smoothed importance sampling. \emph{Journal of Machine Learning Research},
25(72):1-58.
\href{https://jmlr.org/papers/v25/19-556.html}{PDF}

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018) Using
stacking to average Bayesian predictive distributions.
\emph{Bayesian Analysis}, advance publication,  doi:10.1214/17-BA1091.
(\href{https://projecteuclid.org/euclid.ba/1516093227}{online}).

Magnusson, M., Riis Andersen, M., Jonasson, J. and Vehtari, A. (2019).
Leave-One-Out Cross-Validation for Large Data.
In \emph{Thirty-sixth International Conference on Machine Learning},
PMLR 97:4244-4253.

Magnusson, M., Riis Andersen, M., Jonasson, J. and Vehtari, A. (2020).
Leave-One-Out Cross-Validation for Model Comparison in Large Data.
In \emph{Proceedings of the 23rd International Conference on Artificial
Intelligence and Statistics (AISTATS)}, PMLR 108:341-351.

Epifani, I., MacEachern, S. N., and Peruggia, M. (2008). Case-deletion
importance sampling estimators: Central limit theorems and related results.
\emph{Electronic Journal of Statistics} \strong{2}, 774-806.

Gelfand, A. E. (1996). Model determination using sampling-based methods. In
\emph{Markov Chain Monte Carlo in Practice}, ed. W. R. Gilks, S. Richardson,
D. J. Spiegelhalter, 145-162. London: Chapman and Hall.

Gelfand, A. E., Dey, D. K., and Chang, H. (1992). Model determination using
predictive distributions with implementation via sampling-based methods. In
\emph{Bayesian Statistics 4}, ed. J. M. Bernardo, J. O. Berger, A. P. Dawid,
and A. F. M. Smith, 147-167. Oxford University Press.

Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive
information criteria for Bayesian models. \emph{Statistics and Computing}
\strong{24}, 997-1016.

Ionides, E. L. (2008). Truncated importance sampling. \emph{Journal of
Computational and Graphical Statistics} \strong{17}, 295-311.

Koopman, S. J., Shephard, N., and Creal, D. (2009). Testing the assumptions
behind importance sampling. \emph{Journal of Econometrics} \strong{149}, 2-11.

Peruggia, M. (1997). On the variability of case-deletion importance sampling
weights in the Bayesian linear model. \emph{Journal of the American
Statistical Association} \strong{92}, 199-207.

Stan Development Team (2017). The Stan C++ Library, Version 2.17.0.
\url{https://mc-stan.org}.

Stan Development Team (2018). RStan: the R interface to Stan, Version 2.17.3.
\url{https://mc-stan.org}.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
widely application information criterion in singular learning theory.
\emph{Journal of Machine Learning Research} \strong{11}, 3571-3594.

Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation method
for the generalized Pareto distribution. \emph{Technometrics} \strong{51},
316-325.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://mc-stan.org/loo/}
  \item \url{https://discourse.mc-stan.org}
  \item Report bugs at \url{https://github.com/stan-dev/loo/issues}
}

}
\author{
\strong{Maintainer}: Jonah Gabry \email{jsg2201@columbia.edu}

Authors:
\itemize{
  \item Aki Vehtari \email{Aki.Vehtari@aalto.fi}
  \item Måns Magnusson
  \item Yuling Yao
  \item Paul-Christian Bürkner
  \item Topi Paananen
  \item Andrew Gelman
}

Other contributors:
\itemize{
  \item Ben Goodrich [contributor]
  \item Juho Piironen [contributor]
  \item Bruno Nicenboim [contributor]
  \item Leevi Lindgren [contributor]
}

}
