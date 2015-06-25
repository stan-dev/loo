[![Travis-CI Build Status](https://travis-ci.org/jgabry/loo.svg?branch=develop)](https://travis-ci.org/jgabry/loo)

**loo** R package

Efficient leave-one-out cross-validation and WAIC for fitted Bayesian models

### About 

Leave-one-out cross-validation (LOO) and the widely applicable information
criterion (WAIC) are methods for estimating pointwise out-of-sample
prediction accuracy from a fitted Bayesian model using the log-likelihood
evaluated at the posterior simulations of the parameter values. LOO and WAIC
have various advantages over simpler estimates of predictive error such as
AIC and DIC but are less used in practice because they involve additional
computational steps. 

This package implements the fast and stable computations
for LOO and WAIC laid out in our forthcoming paper *Efficient leave-one-out cross-validation and WAIC for evaluating fitted Bayesian models*. From existing posterior simulation draws, we compute LOO using very good importance
sampling (VGIS), a new procedure for regularizing importance weights. As a
byproduct of our calculations, we also obtain approximate standard errors for
estimated predictive errors and for comparing of predictive errors between
two models.

### Authors
Aki Vehtari, Andrew Gelman, Jonah Gabry
