#' Efficient implementation of leave-one-out cross-validation and WAIC for
#' evaluating fitted Bayesian models
#'
#' @description Leave-one-out cross-validation and WAIC
#'
#' Leave-one-out cross-validation (LOO) and the widely applicable information
#' criterion (WAIC) are methods for estimating pointwise out-of-sample
#' prediction accuracy from a fitted Bayesian model using the log-likelihood
#' evaluated at the posterior simulations of the parameter values. LOO and WAIC
#' have various advantages over simpler estimates of predictive error such as
#' AIC and DIC but are less used in practice because they involve additional
#' computational steps. Here we lay out fast and stable computations for LOO and
#' WAIC that can be performed using existing simulation draws. We compute LOO
#' using very good importance sampling (VGIS), a new procedure for regularizing
#' importance weights. As a byproduct of our calculations, we also obtain
#' approximate standard errors for estimated predictive errors and for comparing
#' of predictive errors between two models.
#'
#' @docType package
#' @name loo
#'
NULL
