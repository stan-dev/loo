[![Travis-CI Build Status](https://travis-ci.org/jgabry/loo.svg?branch=master)](https://travis-ci.org/jgabry/loo)
[![codecov.io](https://codecov.io/github/jgabry/loo/coverage.svg?branch=master)](https://codecov.io/github/jgabry/loo?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/loo?color=blue)](http://cran.r-project.org/web/packages/loo)
[![RStudio_CRAN_mirror_downloads_badge](http://cranlogs.r-pkg.org/badges/grand-total/loo?color=blue)](http://cran.r-project.org/web/packages/loo)

### **loo** R package

Efficient leave-one-out cross-validation and WAIC for fitted Bayesian models

### Install

* **CRAN** `install.packages("loo")`

* **GitHub** `devtools::install_github("jgabry/loo")` 

### About 

Leave-one-out cross-validation (LOO) and the widely applicable information
criterion (WAIC) are methods for estimating pointwise out-of-sample
prediction accuracy from a fitted Bayesian model using the log-likelihood
evaluated at the posterior simulations of the parameter values. LOO and WAIC
have various advantages over simpler estimates of predictive error such as
AIC and DIC but are less used in practice because they involve additional
computational steps. 

This package implements the fast and stable computations for LOO and WAIC from
[*Efficient leave-one-out cross-validation and WAIC for evaluating fitted Bayesian models*](http://arxiv.org/abs/1507.04544). 
From existing posterior simulation draws, we compute LOO using Pareto smoothed
importance sampling (PSIS), a new procedure for regularizing importance weights.
As a byproduct of our calculations, we also obtain approximate standard errors
for estimated predictive errors and for comparing predictive errors between two
models.

### Authors
Aki Vehtari, Andrew Gelman, Jonah Gabry

### Python and Matlab/Octave Code
Corresponding Python and Matlab/Octave code can be found at https://github.com/avehtari/PSIS

