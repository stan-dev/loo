#' LOO package glossary
#'
#' @name loo-glossary
#'
#' @template loo-and-psis-references
#' @template bayesvis-reference
#'
#' @description
#'   The pages provides definitions to key terms. Also see the
#'   [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html) on
#'   the __loo__ website for answers to frequently asked questions.
#'
#'   Note: VGG2017 refers to Vehtari, Gelman, and Gabry (2017). See
#'   **References**, below.
#'
#' @section ELPD and `elpd_loo`:
#'
#' The ELPD is the theoretical expected log pointwise predictive density for a new
#' dataset (Eq 1 in VGG2017), which can be estimated, e.g., using
#' cross-validation. `elpd_loo` is the Bayesian LOO estimate of the
#' expected log pointwise predictive density (Eq 4 in VGG2017) and
#' is a sum of N individual pointwise log predictive densities. Probability
#' densities can be smaller or larger than 1, and thus log predictive densities
#' can be negative or positive. For simplicity the ELPD acronym is used also for
#' expected log pointwise predictive probabilities for discrete models.
#' Probabilities are always equal or less than 1, and thus log predictive
#' probabilities are 0 or negative.
#'
#' @section Standard error of `elpd_loo`:
#'
#' As `elpd_loo` is defined as the sum of N independent components (Eq 4 in
#' VGG2017), we can compute the standard error by using the standard deviation
#' of the N components and multiplying by `sqrt(N)` (Eq 23 in VGG2017).
#' This standard error is a coarse description of our uncertainty about the
#' predictive performance for unknown future data. When N is small or there is
#' severe model misspecification, the current SE estimate is overoptimistic and
#' the actual SE can even be twice as large. Even for moderate N, when the SE
#' estimate is an accurate estimate for the scale, it ignores the skewness. When
#' making model comparisons, the SE of the component-wise (pairwise) differences
#' should be used instead (see the `se_diff` section below and Eq 24 in
#' VGG2017).
#'
#' @section Monte Carlo SE of elpd_loo:
#'
#' The Monte Carlo standard error is the estimate for the computational accuracy
#' of MCMC and importance sampling used to compute `elpd_loo`. Usually this
#' is negligible compared to the standard describing the uncertainty due to
#' finite number of observations (Eq 23 in VGG2017).
#'
#' @section `p_loo` (effective number of parameters):
#'
#' `p_loo` is the difference between `elpd_loo` and the non-cross-validated
#' log posterior predictive density. It describes how much more difficult it
#' is to predict future data than the observed data. Asymptotically under
#' certain regularity conditions, `p_loo` can be interpreted as the
#' *effective number of parameters*. In well behaving cases `p_loo < N` and
#' `p_loo < p`, where `p` is the total number of parameters in the
#' model. `p_loo > N`  or `p_loo > p` indicates that the model has very
#' weak predictive capability and may indicate a severe model misspecification.
#' See below for more on interpreting `p_loo` when there are warnings
#' about high Pareto k diagnostic values.
#'
#' @section Pareto k estimates:
#'
#' The Pareto `k` estimate is a diagnostic for Pareto smoothed importance
#' sampling (PSIS), which is used to compute components of `elpd_loo`. In
#' importance-sampling LOO (the full posterior distribution is used as the
#' proposal distribution). The Pareto k diagnostic estimates how far an
#' individual leave-one-out distribution is from the full distribution. If
#' leaving out an observation changes the posterior too much then importance
#' sampling is not able to give reliable estimate. If `k<0.5`, then the
#' corresponding component of `elpd_loo` is estimated with high accuracy.
#' If `0.5<k<0.7` the accuracy is lower, but still ok. If `k>0.7`,
#' then importance sampling is not able to provide useful estimate for that
#' component/observation. Pareto k is also useful as a measure of influence of
#' an observation. Highly influential observations have high k values. Very high
#' k values often indicate model misspecification, outliers or mistakes in data
#' processing. See Section 6 of Gabry et al. (2019) for an example.
#'
#' \subsection{Interpreting `p_loo` when Pareto `k` is large}{
#' If `k > 0.7` then we can also look at the `p_loo` estimate for
#' some additional information about the problem:
#'
#' \itemize{
#' \item If `p_loo << p` (the total number of parameters in the model),
#' then the model is likely to be misspecified. Posterior predictive checks
#' (PPCs) are then likely to also detect the problem. Try using an overdispersed
#' model, or add more structural information (nonlinearity, mixture model,
#' etc.).
#'
#' \item If `p_loo < p` and the number of parameters `p` is relatively
#' large compared to the number of observations (e.g., `p>N/5`), it is
#' likely that the model is so flexible or the population prior so weak that it’s
#' difficult to predict the left out observation (even for the true model).
#' This happens, for example, in the simulated 8 schools (in VGG2017), random
#' effect models with a few observations per random effect, and Gaussian
#' processes and spatial models with short correlation lengths.
#'
#' \item If `p_loo > p`, then the model is likely to be badly misspecified.
#' If the number of parameters `p<<N`, then PPCs are also likely to detect the
#' problem. See the case study at
#' <https://avehtari.github.io/modelselection/roaches.html> for an example.
#' If `p` is relatively large compared to the number of
#' observations, say `p>N/5` (more accurately we should count number of
#' observations influencing each parameter as in hierarchical models some groups
#' may have few observations and other groups many), it is possible that PPCs won't
#' detect the problem.
#' }
#' }
#'
#' @section elpd_diff:
#' `elpd_diff` is the difference in `elpd_loo` for two models. If more
#' than two models are compared, the difference is computed relative to the
#' model with highest `elpd_loo`.
#'
#' @section se_diff:
#'
#' The standard error of component-wise differences of elpd_loo (Eq 24 in
#' VGG2017) between two models. This SE is *smaller* than the SE for
#' individual models due to correlation (i.e., if some observations are easier
#' and some more difficult to predict for all models).
#'
NULL
