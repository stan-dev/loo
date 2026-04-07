#' LOO package glossary
#'
#' @name loo-glossary
#'
#' @template loo-and-psis-references
#' @template loo-uncertainty-reference
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
#' VGG2017). Sivula et al. (2022) discuss the conditions when the normal
#' approximation used for SE and `se_diff` is good.
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
#' The Pareto \eqn{k} estimate is a diagnostic for Pareto smoothed importance
#' sampling (PSIS), which is used to compute components of `elpd_loo`. In
#' importance-sampling LOO the full posterior distribution is used as the
#' proposal distribution. The Pareto k diagnostic estimates how far an
#' individual leave-one-out distribution is from the full distribution. If
#' leaving out an observation changes the posterior too much then importance
#' sampling is not able to give a reliable estimate. Pareto smoothing stabilizes
#' importance sampling and guarantees a finite variance estimate at the
#' cost of some bias.
#'
#' The diagnostic threshold for Pareto \eqn{k} depends on sample size
#' \eqn{S} (sample size dependent threshold was introduced by Vehtari
#' et al., 2024, and before that fixed thresholds of 0.5 and 0.7 were
#' recommended). For simplicity, `loo` package uses the nominal sample
#' size \eqn{S}  when computing the sample size specific
#' threshold. This provides an optimistic threshold if the effective
#' sample size is less than 2200, but even then if ESS/S > 1/2 the difference
#' is usually negligible. Thinning of MCMC draws can be used to improve
#' the ratio ESS/S.
#'
#' * If \eqn{k < \min(1 - 1 / \log_{10}(S), 0.7)}, where \eqn{S} is the
#'   sample size, the PSIS estimate and the corresponding Monte
#'   Carlo standard error estimate are reliable.
#'
#' * If \eqn{1 - 1 / \log_{10}(S) <= k < 0.7}, the PSIS estimate and the
#'   corresponding Monte Carlo standard error estimate are not
#'   reliable, but increasing the (effective) sample size \eqn{S} above
#'   2200 may help (this will increase the sample size specific
#'   threshold \eqn{(1 - 1 / \log_{10}(2200) > 0.7} and then the bias specific
#'   threshold 0.7 dominates).
#'
#' * If \eqn{0.7 <= k < 1}, the PSIS estimate and the corresponding Monte
#'   Carlo standard error have large bias and are not reliable. Increasing
#'   the sample size may reduce the variability in the \eqn{k} estimate, which
#'   may also result in a lower \eqn{k} estimate.
#'
#' * If \eqn{k \geq 1}{k >= 1}, the target distribution is estimated to
#'   have non-finite mean. The PSIS estimate and the corresponding Monte
#'   Carlo standard error are not well defined. Increasing the sample size
#'   may reduce the variability in \eqn{k} estimate, which may also result in
#'   a lower \eqn{k} estimate.
#'
#' Pareto \eqn{k} is also useful as a measure of influence of an
#' observation.  Highly influential observations have high \eqn{k}
#' values. Very high \eqn{k} values often indicate model
#' misspecification, outliers or mistakes in data processing. See
#' Section 6 of Gabry et al. (2019) for an example.
#'
#' \subsection{Interpreting `p_loo` when Pareto `k` is large}{
#' If \eqn{k > 0.7} then we can also look at
#' the `p_loo` estimate for some additional information about the problem:
#'
#' * If `p_loo << p` (the total number of parameters in the model),
#' then the model is likely to be misspecified. Posterior predictive checks
#' (PPCs) are then likely to also detect the problem. Try using an overdispersed
#' model, or add more structural information (nonlinearity, mixture model,
#' etc.).
#'
#' * If `p_loo < p` and the number of parameters `p` is relatively
#' large compared to the number of observations (e.g., `p>N/5`), it is
#' likely that the model is so flexible or the population prior so weak that it’s
#' difficult to predict the left out observation (even for the true model).
#' This happens, for example, in the simulated 8 schools (in VGG2017), random
#' effect models with a few observations per random effect, and Gaussian
#' processes and spatial models with short correlation lengths.
#'
#' * If `p_loo > p`, then the model is likely to be badly misspecified.
#' If the number of parameters `p<<N`, then PPCs are also likely to detect the
#' problem. See the case study at
#' <https://avehtari.github.io/modelselection/roaches.html> for an example.
#' If `p` is relatively large compared to the number of
#' observations, say `p>N/5` (more accurately we should count number of
#' observations influencing each parameter as in hierarchical models some groups
#' may have few observations and other groups many), it is possible that PPCs won't
#' detect the problem.
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
#' @section `p_worse` (probability of worse predictive performance):
#'
#' `p_worse` is the estimated probability that a model has worse predictive
#' performance than the best-ranked model in the comparison, based on the normal
#' approximation to the uncertainty in `elpd_diff`. It is computed as
#'
#'     p_worse = pnorm(0, elpd_diff, se_diff).
#'
#' The best-ranked model (the first row in the `loo_compare()` output, where
#' `elpd_diff = 0`) always receives `NA`, since the comparison is defined
#' relative to that model.
#'
#' Because models are ordered by `elpd_loo` before computing `p_worse`, all
#' reported values are at least 0.5 by construction. A value close to 0.5
#' indicates that the models are nearly indistinguishable in predictive
#' performance and that the ranking could easily be reversed with different
#' data. A value close to 1 indicates that the lower-ranked model is almost
#' certainly worse. `p_worse` inherits all the limitations of `se_diff` and the
#' normal approximation on which it is based. In particular, when `se_diff` is
#' underestimated, `p_worse` will be estimated too close to 1, making a model
#' appear more clearly worse than the data actually support. Conversely, when
#' `elpd_diff` is biased due to an unreliable LOO approximation, `p_worse` can
#' point in the wrong direction entirely. When any of these conditions are
#' present, `diag_diff` or `diag_elpd` will be flagged in the `loo_compare()`
#' output. See those sections below for further guidance.
#'
#' @section `diag_diff` (pairwise comparison diagnostics):
#'
#' `diag_diff` is a diagnostic column in the `loo_compare()` output for each
#' model comparison against the current reference model. It flags conditions
#' under which the normal approximation behind `se_diff` and `p_worse` is likely
#' to be poorly calibrated. The column contains a short label when a condition
#' is detected, and is empty otherwise.
#'
#' The column `diag_diff` currently flags two problems:
#'
#' ### `N < 100`
#'
#' When the number of observations is small, we may assume `se_diff` to be
#' underestimated. As a rough heuristic one can multiply `se_diff` by 2 to
#' make a more conservative estimate.
#'
#' ### `|elpd_diff| < 4`
#'
#' When `|elpd_diff|` is below 4, the models have very similar predictive
#' performance. In this setting, Sivula et al. (2025) show that skewness in
#' the error distribution can make the normal approximation for `se_diff`
#' and `p_worse` miscalibrated, even for large N. In practice, this usually
#' supports treating the models as predictively similar.
#'
#' ### Relation between `N < 100` and `|elpd_diff| < 4`
#'
#' The conditions flagged by `diag_diff` are not independent: they tend to
#' co-occur, and when they do, some flags carry more information than others.
#' `loo_compare()` therefore follows a priority hierarchy and shows only the
#' most critical flag in the table output.
#'
#' The hierarchy is as follows:
#'
#' * **`N < 100` takes highest priority.** A small sample size undermines the
#' reliability of `se_diff` by underestimating uncertainty. Because of this,
#' even if `|elpd_diff| < 4` is also true for a comparison, the table will only
#' show `N < 100`. The small sample size renders the `|elpd_diff| < 4`
#' diagnostic less meaningful.
#'
#' * **`|elpd_diff| < 4` takes second priority.** When N >= 100 and the
#'   difference is small, the normal approximation is miscalibrated due to the
#'   skewness of the error distribution (Sivula et al., 2025). In this
#'   situation, `se_diff` exists and is not heavily biased in scale, but the
#'   shape of the approximation is wrong, making `p_worse` unreliable.
#'
#' @section `diag_elpd`:
#'
#' `diag_elpd` is a diagnostic column in the `loo_compare()` output that flags
#' when the PSIS-LOO approximation for an individual model is unreliable. Unlike
#' `diag_diff`, which concerns the *comparison* between models, `diag_elpd`
#' concerns the quality of the `elpd_loo` estimate for each model individually.
#' It contains a short text label when a problem is detected, and is empty
#' otherwise.
#'
#' ### `K k_psis > t` (K observations with Pareto-k values > t)
#'
#' This label indicates that K observations for this model have Pareto-k values
#' above the PSIS reliability threshold `t` used by `loo` for that fit. The
#' threshold is sample-size dependent, and in many practical cases close to
#' 0.7. When this flag appears, the PSIS approximation can be unreliable for
#' those observations, and the resulting `elpd_loo` may be biased. Because
#' `elpd_diff` is a direct difference of two models' `elpd_loo` values, bias in
#' either model's estimate propagates directly into `elpd_diff` and `p_worse`.
#' This is qualitatively different from the calibration issues flagged by
#' `diag_diff`: here the estimate itself may be wrong, not just uncertain.
#'
#' See for further information on Pareto-k values the "Pareto k estimates"
#' section.
NULL
