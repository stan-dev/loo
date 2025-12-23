# Package index

## Package description, glossary, and included data sets

- [`loo-package`](https://mc-stan.org/loo/reference/loo-package.md) :
  Efficient LOO-CV and WAIC for Bayesian models
- [`loo-glossary`](https://mc-stan.org/loo/reference/loo-glossary.md) :
  LOO package glossary
- [`loo-datasets`](https://mc-stan.org/loo/reference/loo-datasets.md)
  [`Kline`](https://mc-stan.org/loo/reference/loo-datasets.md)
  [`milk`](https://mc-stan.org/loo/reference/loo-datasets.md)
  [`voice`](https://mc-stan.org/loo/reference/loo-datasets.md)
  [`voice_loo`](https://mc-stan.org/loo/reference/loo-datasets.md) :
  Datasets for loo examples and vignettes

## Approximate LOO-CV

Approximate LOO-CV, Pareto smoothed importance sampling (PSIS), and
diagnostics.

- [`loo()`](https://mc-stan.org/loo/reference/loo.md)
  [`loo_i()`](https://mc-stan.org/loo/reference/loo.md)
  [`is.loo()`](https://mc-stan.org/loo/reference/loo.md)
  [`is.psis_loo()`](https://mc-stan.org/loo/reference/loo.md) :
  Efficient approximate leave-one-out cross-validation (LOO)
- [`loo_subsample()`](https://mc-stan.org/loo/reference/loo_subsample.md)
  : Efficient approximate leave-one-out cross-validation (LOO) using
  subsampling, so that less costly and more approximate computation is
  made for all LOO-fold, and more costly and accurate computations are
  made only for m\<N LOO-folds.
- [`loo_approximate_posterior()`](https://mc-stan.org/loo/reference/loo_approximate_posterior.md)
  : Efficient approximate leave-one-out cross-validation (LOO) for
  posterior approximations
- [`loo_moment_match()`](https://mc-stan.org/loo/reference/loo_moment_match.md)
  : Moment matching for efficient approximate leave-one-out
  cross-validation (LOO)
- [`loo_moment_match_split()`](https://mc-stan.org/loo/reference/loo_moment_match_split.md)
  : Split moment matching for efficient approximate leave-one-out
  cross-validation (LOO)
- [`E_loo()`](https://mc-stan.org/loo/reference/E_loo.md) : Compute
  weighted expectations
- [`psis()`](https://mc-stan.org/loo/reference/psis.md)
  [`is.psis()`](https://mc-stan.org/loo/reference/psis.md)
  [`is.sis()`](https://mc-stan.org/loo/reference/psis.md)
  [`is.tis()`](https://mc-stan.org/loo/reference/psis.md) : Pareto
  smoothed importance sampling (PSIS)
- [`ap_psis()`](https://mc-stan.org/loo/reference/ap_psis.md) : Pareto
  smoothed importance sampling (PSIS) using approximate posteriors
- [`tis()`](https://mc-stan.org/loo/reference/tis.md) : Truncated
  importance sampling (TIS)
- [`sis()`](https://mc-stan.org/loo/reference/sis.md) : Standard
  importance sampling (SIS)
- [`importance_sampling()`](https://mc-stan.org/loo/reference/importance_sampling.md)
  : A parent class for different importance sampling methods.
- [`weights(`*`<importance_sampling>`*`)`](https://mc-stan.org/loo/reference/weights.importance_sampling.md)
  : Extract importance sampling weights
- [`pareto_k_table()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  [`pareto_k_ids()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  [`pareto_k_values()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  [`pareto_k_influence_values()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  [`psis_n_eff_values()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  [`mcse_loo()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  [`plot(`*`<psis_loo>`*`)`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  [`plot(`*`<psis>`*`)`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
  : Diagnostics for Pareto smoothed importance sampling (PSIS)

## Model comparison weighting/averaging

Functions for comparing models and computing model weights via stacking
of predictive distributions or pseudo-BMA weighting.

- [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.md)
  [`print(`*`<compare.loo>`*`)`](https://mc-stan.org/loo/reference/loo_compare.md)
  [`print(`*`<compare.loo_ss>`*`)`](https://mc-stan.org/loo/reference/loo_compare.md)
  : Model comparison
- [`loo_model_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.md)
  [`stacking_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.md)
  [`pseudobma_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.md)
  : Model averaging/weighting via stacking or pseudo-BMA weighting

## Helper functions for K-fold CV

- [`kfold_split_random()`](https://mc-stan.org/loo/reference/kfold-helpers.md)
  [`kfold_split_stratified()`](https://mc-stan.org/loo/reference/kfold-helpers.md)
  [`kfold_split_grouped()`](https://mc-stan.org/loo/reference/kfold-helpers.md)
  : Helper functions for K-fold cross-validation
- [`kfold()`](https://mc-stan.org/loo/reference/kfold-generic.md)
  [`is.kfold()`](https://mc-stan.org/loo/reference/kfold-generic.md) :
  Generic function for K-fold cross-validation for developers
- [`elpd()`](https://mc-stan.org/loo/reference/elpd.md) : Generic
  (expected) log-predictive density

## Other functions

- [`loo_predictive_metric()`](https://mc-stan.org/loo/reference/loo_predictive_metric.md)
  : Estimate leave-one-out predictive performance..

- [`crps()`](https://mc-stan.org/loo/reference/crps.md)
  [`scrps()`](https://mc-stan.org/loo/reference/crps.md)
  [`loo_crps()`](https://mc-stan.org/loo/reference/crps.md)
  [`loo_scrps()`](https://mc-stan.org/loo/reference/crps.md) :
  Continuously ranked probability score

- [`elpd()`](https://mc-stan.org/loo/reference/elpd.md) : Generic
  (expected) log-predictive density

- [`waic()`](https://mc-stan.org/loo/reference/waic.md)
  [`is.waic()`](https://mc-stan.org/loo/reference/waic.md) : Widely
  applicable information criterion (WAIC)

- [`extract_log_lik()`](https://mc-stan.org/loo/reference/extract_log_lik.md)
  : Extract pointwise log-likelihood from a Stan model

- [`pointwise()`](https://mc-stan.org/loo/reference/pointwise.md) :
  Convenience function for extracting pointwise estimates

- [`relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.md)
  : Convenience function for computing relative efficiencies

- [`gpdfit()`](https://mc-stan.org/loo/reference/gpdfit.md) : Estimate
  parameters of the Generalized Pareto distribution

- [`example_loglik_array()`](https://mc-stan.org/loo/reference/example_loglik_array.md)
  [`example_loglik_matrix()`](https://mc-stan.org/loo/reference/example_loglik_array.md)
  : Objects to use in examples and tests

- [`print(`*`<loo>`*`)`](https://mc-stan.org/loo/reference/print.loo.md)
  [`print(`*`<waic>`*`)`](https://mc-stan.org/loo/reference/print.loo.md)
  [`print(`*`<psis_loo>`*`)`](https://mc-stan.org/loo/reference/print.loo.md)
  [`print(`*`<importance_sampling_loo>`*`)`](https://mc-stan.org/loo/reference/print.loo.md)
  [`print(`*`<psis_loo_ap>`*`)`](https://mc-stan.org/loo/reference/print.loo.md)
  [`print(`*`<psis>`*`)`](https://mc-stan.org/loo/reference/print.loo.md)
  [`print(`*`<importance_sampling>`*`)`](https://mc-stan.org/loo/reference/print.loo.md)
  : Print methods

- [`nobs(`*`<psis_loo_ss>`*`)`](https://mc-stan.org/loo/reference/nobs.psis_loo_ss.md)
  :

  The number of observations in a `psis_loo_ss` object.

- [`obs_idx()`](https://mc-stan.org/loo/reference/obs_idx.md) : Get
  observation indices used in subsampling

- [`update(`*`<psis_loo_ss>`*`)`](https://mc-stan.org/loo/reference/update.psis_loo_ss.md)
  :

  Update `psis_loo_ss` objects

## Deprecated functions

- [`compare()`](https://mc-stan.org/loo/reference/compare.md) : Model
  comparison (deprecated, old version)
- [`psislw()`](https://mc-stan.org/loo/reference/psislw.md) : Pareto
  smoothed importance sampling (deprecated, old version)
