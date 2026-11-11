# Package index

## Package description, glossary, and included data sets

- [`loo-package`](https://mc-stan.org/loo/pr/363/reference/loo-package.md)
  : Efficient LOO-CV and WAIC for Bayesian models
- [`loo-glossary`](https://mc-stan.org/loo/pr/363/reference/loo-glossary.md)
  : LOO package glossary
- [`loo-datasets`](https://mc-stan.org/loo/pr/363/reference/loo-datasets.md)
  [`Kline`](https://mc-stan.org/loo/pr/363/reference/loo-datasets.md)
  [`milk`](https://mc-stan.org/loo/pr/363/reference/loo-datasets.md)
  [`voice`](https://mc-stan.org/loo/pr/363/reference/loo-datasets.md)
  [`voice_loo`](https://mc-stan.org/loo/pr/363/reference/loo-datasets.md)
  : Datasets for loo examples and vignettes

## Approximate LOO-CV

Approximate LOO-CV, Pareto smoothed importance sampling (PSIS), and
diagnostics.

- [`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md)
  [`loo_i()`](https://mc-stan.org/loo/pr/363/reference/loo.md)
  [`is.loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md)
  [`is.psis_loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md) :
  Efficient approximate leave-one-out cross-validation (LOO)
- [`loo_subsample()`](https://mc-stan.org/loo/pr/363/reference/loo_subsample.md)
  : Efficient approximate LOO-CV using subsampling
- [`loo_approximate_posterior()`](https://mc-stan.org/loo/pr/363/reference/loo_approximate_posterior.md)
  : Efficient approximate leave-one-out cross-validation (LOO) for
  posterior approximations
- [`loo_moment_match()`](https://mc-stan.org/loo/pr/363/reference/loo_moment_match.md)
  : Moment matching for efficient approximate leave-one-out
  cross-validation (LOO)
- [`loo_moment_match_split()`](https://mc-stan.org/loo/pr/363/reference/loo_moment_match_split.md)
  : Split moment matching for efficient approximate leave-one-out
  cross-validation (LOO)
- [`E_loo()`](https://mc-stan.org/loo/pr/363/reference/E_loo.md) :
  Compute weighted expectations
- [`psis()`](https://mc-stan.org/loo/pr/363/reference/psis.md)
  [`is.psis()`](https://mc-stan.org/loo/pr/363/reference/psis.md)
  [`is.sis()`](https://mc-stan.org/loo/pr/363/reference/psis.md)
  [`is.tis()`](https://mc-stan.org/loo/pr/363/reference/psis.md) :
  Pareto smoothed importance sampling (PSIS)
- [`ap_psis()`](https://mc-stan.org/loo/pr/363/reference/ap_psis.md) :
  Pareto smoothed importance sampling (PSIS) using approximate
  posteriors
- [`tis()`](https://mc-stan.org/loo/pr/363/reference/tis.md) : Truncated
  importance sampling (TIS)
- [`sis()`](https://mc-stan.org/loo/pr/363/reference/sis.md) : Standard
  importance sampling (SIS)
- [`importance_sampling()`](https://mc-stan.org/loo/pr/363/reference/importance_sampling.md)
  : A parent class for different importance sampling methods.
- [`weights(`*`<importance_sampling>`*`)`](https://mc-stan.org/loo/pr/363/reference/weights.importance_sampling.md)
  : Extract importance sampling weights
- [`pareto_k_table()`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  [`pareto_k_ids()`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  [`pareto_k_values()`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  [`pareto_k_influence_values()`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  [`psis_n_eff_values()`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  [`mcse_loo()`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  [`plot(`*`<psis_loo>`*`)`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  [`plot(`*`<psis>`*`)`](https://mc-stan.org/loo/pr/363/reference/pareto-k-diagnostic.md)
  : Diagnostics for Pareto smoothed importance sampling (PSIS)

## Model comparison weighting/averaging

Functions for comparing models and computing model weights via stacking
of predictive distributions or pseudo-BMA weighting.

- [`print(`*`<compare.loo>`*`)`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
  [`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
  : Model comparison
- [`loo_compare()`](https://mc-stan.org/loo/pr/363/reference/loo_compare.md)
  : Model comparison (deprecated)
- [`loo_model_weights()`](https://mc-stan.org/loo/pr/363/reference/loo_model_weights.md)
  [`stacking_weights()`](https://mc-stan.org/loo/pr/363/reference/loo_model_weights.md)
  [`pseudobma_weights()`](https://mc-stan.org/loo/pr/363/reference/loo_model_weights.md)
  : Model averaging/weighting via stacking or pseudo-BMA weighting

## Helper functions for K-fold CV

- [`kfold_split_random()`](https://mc-stan.org/loo/pr/363/reference/kfold-helpers.md)
  [`kfold_split_stratified()`](https://mc-stan.org/loo/pr/363/reference/kfold-helpers.md)
  [`kfold_split_grouped()`](https://mc-stan.org/loo/pr/363/reference/kfold-helpers.md)
  : Helper functions for K-fold cross-validation
- [`kfold()`](https://mc-stan.org/loo/pr/363/reference/kfold-generic.md)
  [`is.kfold()`](https://mc-stan.org/loo/pr/363/reference/kfold-generic.md)
  : Generic function for K-fold cross-validation for developers

## Predictive performance measures

Unified API for predictive performance measures (density scores, scoring
rules, and point-prediction metrics) with in-sample, LOO, K-fold, and
holdout evaluation modes.

- [`insample_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/insample_pred_measure.md)
  : In-sample predictive performance measures

- [`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)
  : PSIS-LOO predictive performance measures

- [`kfold_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/kfold_pred_measure.md)
  : K-fold cross-validation predictive performance measures

- [`test_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/test_pred_measure.md)
  : Holdout predictive performance measures

- [`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md)
  : Add predictive performance measures to an existing result

- [`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md)
  : Define a custom predictive measure

- [`supported_measures_list`](https://mc-stan.org/loo/pr/363/reference/supported_measures_list.md)
  : Supported predictive measure names

- [`measure_acc()`](https://mc-stan.org/loo/pr/363/reference/measure_acc.md)
  :

  Classification accuracy (`acc`)

- [`measure_bacc()`](https://mc-stan.org/loo/pr/363/reference/measure_bacc.md)
  :

  Balanced classification accuracy (`bacc`)

- [`measure_brier()`](https://mc-stan.org/loo/pr/363/reference/measure_brier.md)
  :

  Brier score (`brier`)

- [`measure_elpd()`](https://mc-stan.org/loo/pr/363/reference/measure_elpd.md)
  :

  Expected log pointwise predictive density (`elpd`)

- [`measure_ic()`](https://mc-stan.org/loo/pr/363/reference/measure_ic.md)
  :

  Information Criteria (`ic`)

- [`measure_mae()`](https://mc-stan.org/loo/pr/363/reference/measure_mae.md)
  :

  Mean absolute error (`mae`)

- [`measure_mlpd()`](https://mc-stan.org/loo/pr/363/reference/measure_mlpd.md)
  :

  Mean log pointwise predictive density (`mlpd`)

- [`measure_mse()`](https://mc-stan.org/loo/pr/363/reference/measure_mse.md)
  :

  Mean squared error (`mse`)

- [`measure_r2()`](https://mc-stan.org/loo/pr/363/reference/measure_r2.md)
  :

  Predictive R-squared (`r2`)

- [`measure_rmse()`](https://mc-stan.org/loo/pr/363/reference/measure_rmse.md)
  :

  Root mean squared error (`rmse`)

- [`measure_rps()`](https://mc-stan.org/loo/pr/363/reference/measure_rps.md)
  : Ranked Probability Score (RPS, SRPS, CRPS, SCRPS)

- [`measure_srps()`](https://mc-stan.org/loo/pr/363/reference/measure_srps.md)
  : Scaled Ranked Probability Score (SRPS, SCRPS)

- [`ptw_log_pred_density()`](https://mc-stan.org/loo/pr/363/reference/ptw_log_pred_density.md)
  :

  Pointwise log predictive density (`lppd_i`)

## Other functions

- [`waic()`](https://mc-stan.org/loo/pr/363/reference/waic.md)
  [`is.waic()`](https://mc-stan.org/loo/pr/363/reference/waic.md) :
  Widely applicable information criterion (WAIC)

- [`extract_log_lik()`](https://mc-stan.org/loo/pr/363/reference/extract_log_lik.md)
  : Extract pointwise log-likelihood from a Stan model

- [`pointwise()`](https://mc-stan.org/loo/pr/363/reference/pointwise.md)
  : Convenience function for extracting pointwise estimates

- [`relative_eff()`](https://mc-stan.org/loo/pr/363/reference/relative_eff.md)
  : Convenience function for computing relative efficiencies

- [`find_model_names()`](https://mc-stan.org/loo/pr/363/reference/find_model_names.md)
  :

  Find the model names associated with `"loo"` objects

- [`gpdfit()`](https://mc-stan.org/loo/pr/363/reference/gpdfit.md) :
  Estimate parameters of the Generalized Pareto distribution

- [`example_loglik_array()`](https://mc-stan.org/loo/pr/363/reference/example_loglik_array.md)
  [`example_loglik_matrix()`](https://mc-stan.org/loo/pr/363/reference/example_loglik_array.md)
  : Objects to use in examples and tests

- [`print(`*`<loo>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  [`print(`*`<waic>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  [`print(`*`<psis_loo>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  [`print(`*`<importance_sampling_loo>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  [`print(`*`<psis_loo_ap>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  [`print(`*`<psis>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  [`print(`*`<importance_sampling>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  [`print(`*`<kfold>`*`)`](https://mc-stan.org/loo/pr/363/reference/print.loo.md)
  : Print methods

- [`nobs(`*`<psis_loo_ss>`*`)`](https://mc-stan.org/loo/pr/363/reference/nobs.psis_loo_ss.md)
  :

  The number of observations in a `psis_loo_ss` object.

- [`obs_idx()`](https://mc-stan.org/loo/pr/363/reference/obs_idx.md) :
  Get observation indices used in subsampling

- [`update(`*`<psis_loo_ss>`*`)`](https://mc-stan.org/loo/pr/363/reference/update.psis_loo_ss.md)
  :

  Update `psis_loo_ss` objects

## Deprecated functions

- [`elpd()`](https://mc-stan.org/loo/pr/363/reference/elpd.md) : Generic
  (expected) log-predictive density (deprecated)
- [`crps()`](https://mc-stan.org/loo/pr/363/reference/crps.md)
  [`scrps()`](https://mc-stan.org/loo/pr/363/reference/crps.md)
  [`loo_crps()`](https://mc-stan.org/loo/pr/363/reference/crps.md)
  [`loo_scrps()`](https://mc-stan.org/loo/pr/363/reference/crps.md) :
  Continuously ranked probability score (deprecated)
- [`loo_predictive_metric()`](https://mc-stan.org/loo/pr/363/reference/loo_predictive_metric.md)
  : Estimate leave-one-out predictive performance (deprecated)
- [`compare()`](https://mc-stan.org/loo/pr/363/reference/compare.md) :
  Model comparison (deprecated, old version)
- [`psislw()`](https://mc-stan.org/loo/pr/363/reference/psislw.md) :
  Pareto smoothed importance sampling (deprecated, old version)
