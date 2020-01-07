[<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo_tm.png" width=100 alt="Stan Logo"/>](https://mc-stan.org)

# loo

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/loo?color=blue)](https://cran.r-project.org/web/packages/loo)
[![RStudio_CRAN_mirror_downloads_badge](https://cranlogs.r-pkg.org/badges/loo?color=blue)](https://cran.r-project.org/web/packages/loo)
[![codecov](https://codecov.io/gh/stan-dev/loo/branch/master/graph/badge.svg)](https://codecov.io/github/stan-dev/loo?branch=master)
[![Travis-CI Build Status](https://travis-ci.org/stan-dev/loo.svg?branch=master)](https://travis-ci.org/stan-dev/loo)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/stan-dev/loo?branch=master&svg=true)](https://ci.appveyor.com/project/jgabry/loo)
<!-- badges: end -->

### Efficient approximate leave-one-out cross-validation for fitted Bayesian models

__loo__ is an R package that allows users to compute efficient approximate
leave-one-out cross-validation for fitted Bayesian models, as well as model
weights that can be used to average predictive distributions.

Leave-one-out cross-validation (LOO-CV, or LOO for short) and the widely
applicable information criterion (WAIC) are methods for estimating pointwise
out-of-sample prediction accuracy from a fitted Bayesian model using the
log-likelihood evaluated at the posterior simulations of the parameter values.
LOO and WAIC have various advantages over simpler estimates of predictive error
such as AIC and DIC but are less used in practice because they involve
additional computational steps.

The __loo__ R package package implements the fast and stable computations 
for approximate LOO-CV and WAIC from

* Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model 
evaluation using leave-one-out cross-validation and WAIC. 
_Statistics and Computing_. 27(5), 1413--1432. 
doi:10.1007/s11222-016-9696-4. [Online](https://link.springer.com/article/10.1007/s11222-016-9696-4), 
[arXiv preprint arXiv:1507.04544](https://arxiv.org/abs/1507.04544).

* Vehtari, A., Gelman, A., and Gabry, J. (2017). Pareto smoothed importance sampling. 
[arXiv preprint arXiv:1507.02646](https://arxiv.org/abs/1507.02646).

From existing posterior simulation draws, we compute approximate LOO-CV using
Pareto smoothed importance sampling (PSIS), a new procedure for regularizing
importance weights. As a byproduct of our calculations, we also obtain
approximate standard errors for estimated predictive errors and for comparing
predictive errors between two models.

We recommend PSIS-LOO-CV instead of WAIC, because PSIS provides useful
diagnostics and effective sample size and Monte Carlo standard error estimates.

As of version `2.0.0`, the __loo__ package also provides methods for using
stacking and other model weighting techiques to average Bayesian predictive
distributions. For details on stacking and model weighting see:

* Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018). Using
stacking to average Bayesian predictive distributions. In Bayesian
Analysis, doi:10.1214/17-BA1091. 
[Online](https://projecteuclid.org/euclid.ba/1516093227),
[arXiv preprint arXiv:1704.02030](https://arxiv.org/abs/1704.02030).


### Resources

* [mc-stan.org/loo](https://mc-stan.org/loo) (online documentation, vignettes)
* [Ask a question](https://discourse.mc-stan.org) (Stan Forums on Discourse)
* [Open an issue](https://github.com/stan-dev/loo/issues) (GitHub issues for bug reports, feature requests)


### Installation

* Install from CRAN:

```r
install.packages("loo")
```

* Install from GitHub (requires __devtools__ package):

```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("stan-dev/loo")
```
We do _not_ recommend setting `build_vignettes=TRUE` when installing from GitHub
because the vignettes take a long time to build and are always available
online at [mc-stan.org/loo/articles/](https://mc-stan.org/loo/articles/).

### Python and Matlab/Octave Code
Corresponding Python and Matlab/Octave code can be found at the
[avehtari/PSIS](https://github.com/avehtari/PSIS) repository.

