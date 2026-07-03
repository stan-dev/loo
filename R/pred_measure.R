#' In-sample predictive performance measures
#'
#' @description
#' Compute predictive performance measures on the same data used to fit the
#' model. This is the simplest entry point when you want density scores
#' (`elpd`) and optional distributional or point-prediction metrics in
#' one call.
#'
#' In-sample `elpd` sums the expected log pointwise predictive density (ELPD)
#' over the training observations. Because the model has already seen these
#' data, in-sample scores are **optimistically biased** for predicting future
#' or otherwise unseen observations. For out-of-sample performance, use
#' [loo_pred_measure()], [kfold_pred_measure()], or [test_pred_measure()]
#' instead; see Vehtari's
#' [Cross-validation FAQ](https://users.aalto.fi/~ave/CV-FAQ.html) (Section 4).
#'
#' @inheritParams pred_measure_params
#'
#' @return
#' An object of class `"insample_pred_measure"` and `"pred_measure"`: a list
#' with:
#' \describe{
#'   \item{`estimates`}{Matrix of summary estimates and standard errors (rows
#'     are measures, columns are `Estimate` and `SE`). The base row `elpd` is
#'     always present when `ylp` is supplied.}
#'   \item{`pointwise`}{Matrix of observation-level contributions (one column
#'     per measure).}
#' }
#'
#' The attribute `source` is `"insample"`. Attribute `dims` gives posterior
#' draws × observations. Use [print()] for a readable summary table.
#'
#' @details
#' **Input requirements by measure.** Supply only the inputs each measure
#' needs:
#'
#' | Measure | `ylp` | `y` | `ypred` | `mupred` |
#' |:---|:---:|:---:|:---:|:---:|
#' | `elpd`, `mlpd`, `ic` | ✓ | | | |
#' | `crps`, `scrps`, `rps`, `srps` | | ✓ | ✓ | |
#' | `acc`, `bacc` | | ✓ | | ✓ |
#' | `mae`, `mse`, `rmse`, `r2` | | ✓ | | ✓ |
#'
#' Base measure `elpd` is always computed when `ylp` is provided. Request
#' `ic`, `mlpd`, or other density scores via `measure`, or supply a custom function;
#' see [supported_measures_list] and the
#' [overview of scores and metrics](https://mc-stan.org/loo/articles/articles-online-only/overview-measures.html)
#' article for definitions and orientation (higher vs lower is better).
#'
#' **Custom measures.** A function passed to `measure` must have attribute
#' `measure_name` and return `estimate`, `se`, and `pointwise`. Only arguments
#' declared in the function signature among `y`, `ypred`, `mupred`, `ylp`, and
#' `log_weights` are supplied automatically.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit <- brms::brm(
#'     Reaction ~ Days, data = lme4::sleepstudy,
#'     refresh = 0, chains = 2, iter = 1000
#'   )
#'   insample_pred_measure(
#'     ylp = brms::log_lik(fit),
#'     y = fit$data$Reaction,
#'     ypred = brms::posterior_predict(fit),
#'     mupred = brms::posterior_epred(fit),
#'     measure = c("rmse", "r2")
#'   )
#' }
#' }
#' \dontrun{
#' # Custom measure (same contract as built-in point-prediction metrics)
#' my_abs_err <- function(y, mupred, log_weights = NULL) {
#'   mu <- colMeans(mupred)
#'   pw <- abs(y - mu)
#'   list(
#'     estimate = mean(pw),
#'     se = sd(pw) / sqrt(length(pw)),
#'     pointwise = pw
#'   )
#' }
#' attr(my_abs_err, "measure_name") <- "my_abs_err"
#' # insample_pred_measure(y = y, mupred = mupred, ylp = ylp, measure = my_abs_err)
#' }
#'
#' @seealso [pred_measure()] to add measures incrementally,
#'   [loo_pred_measure()], [kfold_pred_measure()], [test_pred_measure()],
#'   [supported_measures_list],
#'   [pred-measure workflow article](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html)
#'
#' @export
insample_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measure = NULL,
  group_ids = NULL,
  save_psis = FALSE,
  control = list()
) {
  do_pred_measure(
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    measure = measure, 
    predperf = NULL,
    loo = NULL, 
    kfold = NULL,
    group_ids = group_ids,
    psis_object = NULL,
    save_psis = save_psis,
    source = "insample",
    control = control
  )
}

#' PSIS-LOO predictive performance measures
#'
#' @description
#' Estimate out-of-sample predictive performance with **PSIS-LOO**
#' (Pareto-smoothed importance sampling leave-one-out cross-validation).
#' PSIS-LOO approximates exact LOO-CV without refitting the model once per
#' observation: each held-out point is scored by reweighting the full-data
#' posterior draws.
#'
#' The primary summary is `elpd_loo`, the LOO estimate of expected log
#' pointwise predictive density (ELPD). The base result also includes `p_loo`
#' (effective number of parameters, the difference between in-sample and LOO
#' log predictive density). See [loo::loo()] and the
#' [Cross-validation FAQ](https://users.aalto.fi/~ave/CV-FAQ.html) for
#' interpretation.
#'
#' @inheritParams pred_measure_params
#'
#' @return
#' An object of class `"loo_pred_measure"`, `"pred_measure"`, and (when a
#' `loo` object is supplied) `"loo"`. In addition to `estimates` and
#' `pointwise`, the list may contain:
#' \describe{
#'   \item{`diagnostics`}{PSIS diagnostics, including Pareto \eqn{\hat{k}} in
#'     `diagnostics$pareto_k`. Values above 0.7 suggest unreliable
#'     LOO estimates for those observations.}
#'   \item{`log_weights`}{Normalized log importance weights used for LOO
#'     scoring.}
#'   \item{`psis_object`}{Stored when `save_psis = TRUE`; needed to add
#'     further measures with [pred_measure()] without recomputing weights.}
#' }
#'
#' Measure names carry a `_loo` suffix (e.g. `elpd_loo`, `crps_loo`).
#'
#' @details
#' **Three equivalent input patterns:**
#'
#' \describe{
#'   \item{Precomputed `loo` object}{`loo_pred_measure(loo = loo_fit, ...)`.
#'     Run [loo::loo()] with `save_psis = TRUE`.}
#'   \item{`ylp` + `psis_object`}{Pass both when you have already computed
#'     PSIS weights separately.}
#'   \item{`ylp` only}{PSIS weights are computed internally from `ylp`.}
#' }
#'
#' For distributional and point-prediction measures (`crps`, `r2`, etc.),
#' supply `y`, `ypred`, and/or `mupred` as for [insample_pred_measure()]. When
#' adding measures incrementally, call [pred_measure()] with `predperf` set to
#' an existing result; use `save_psis = TRUE` on the initial call so weights
#' are stored.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit <- brms::brm(
#'     Reaction ~ Days, data = lme4::sleepstudy,
#'     refresh = 0, chains = 2, iter = 1000
#'   )
#'   loo_fit <- loo::loo(fit, save_psis = TRUE)
#'   loo_pred_measure(
#'     loo = loo_fit,
#'     y = fit$data$Reaction,
#'     ypred = brms::posterior_predict(fit),
#'     measure = c("rmse", "r2")
#'   )
#' }
#' }
#'
#' @seealso [insample_pred_measure()], [pred_measure()], [loo::loo()],
#'   [supported_measures_list],
#'   [pred-measure workflow article](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html)
#'
#' @export
loo_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measure = NULL,
  loo = NULL,
  group_ids = NULL,
  psis_object = NULL,
  save_psis = FALSE,
  control = list()
) {
  do_pred_measure(
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_test = NULL,
    measure = measure, 
    predperf = NULL,
    loo = loo, 
    kfold = NULL,
    group_ids = group_ids,
    psis_object = psis_object,
    save_psis = save_psis,
    source = "loo",
    control = control
  )
}

#' K-fold cross-validation predictive performance measures
#'
#' @description
#' Compute predictive performance measures under **k-fold cross-validation**.
#' K-fold CV holds out groups of observations, refits (or reuses stored fits),
#' and scores the held-out folds.
#'
#' Pass a `kfold` object from [brms::kfold()] (with `save_fits = TRUE` when
#' you need posterior predictions on held-out folds). Base density summaries
#' (`elpd_kfold`, `ic_kfold`, `p_kfold`) come from the `kfold` object;
#' additional measures require the same optional inputs as
#' [insample_pred_measure()].
#'
#' @inheritParams pred_measure_params
#'
#' @return
#' An object of class `"kfold_pred_measure"` and `"pred_measure"`, inheriting
#' attributes from the `kfold` object (`K`, `folds`, `fold_type`, etc.). The
#' list contains `estimates` and `pointwise`; measure names carry a `_kfold`
#' suffix (e.g. `elpd_kfold`, `crps_kfold`).
#'
#' @details
#' For distributional measures on held-out folds, obtain posterior predictions
#' with `brms::kfold_predict()` and pass the resulting `yrep` matrices as
#' `ypred` and/or `mupred`. See the sleep-study workflow in
#' [pred-measure workflow article](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit <- brms::brm(
#'     Reaction ~ Days, data = lme4::sleepstudy,
#'     refresh = 0, chains = 2, iter = 1000
#'   )
#'   kf <- brms::kfold(fit, K = 5, save_fits = TRUE)
#'   ypred_kf <- brms::kfold_predict(kf, method = "predict")$yrep
#'   kfold_pred_measure(
#'     y = fit$data$Reaction,
#'     ypred = ypred_kf,
#'     kfold = kf,
#'     measure = "rmse"
#'   )
#' }
#' }
#'
#' @seealso [loo_pred_measure()], [insample_pred_measure()], [pred_measure()],
#'   [brms::kfold()], [supported_measures_list],
#'   [pred-measure workflow article](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html)
#'
#' @export
kfold_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measure = NULL,
  kfold = NULL,
  group_ids = NULL,
  control = list()
) {
  do_pred_measure(
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_test = NULL,
    measure = measure, 
    predperf = NULL,
    loo = NULL, 
    kfold = kfold,
    group_ids = group_ids,
    psis_object = NULL,
    save_psis = FALSE,
    source = "kfold",
    control = control
  )
}

#' Holdout predictive performance measures
#'
#' @description
#' Score predictive performance on **genuinely new (holdout) data** that was
#' not used to fit the model. This mirrors the cross-validation goal of
#' assessing how well a model predicts unseen observations, but with an
#' explicit train/test split rather than LOO or k-fold reweighting.
#'
#' Supply `ylp_test` from log predictive densities evaluated on the holdout set
#' (e.g. `brms::log_lik(fit, newdata = test_data)`). This is required for the
#' base summary `elpd_test`. Optional distributional and point-prediction
#' measures use observed and predicted values on the test set only. Pass training
#' `ylp` only when an additional measure needs log predictive densities from the
#' training fit.
#'
#' @inheritParams pred_measure_params
#'
#' @return
#' An object of class `"test_pred_measure"` and `"pred_measure"` with
#' `estimates` and `pointwise`. Measure names carry a `_test` suffix (e.g.
#' `elpd_test`, `crps_test`). Attribute `dims` reflects the test-set size
#' (from `ylp_test`), not the training data.
#'
#' @details
#' The base summary `elpd_test` is computed from `ylp_test` on the holdout
#' observations only.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   data <- lme4::sleepstudy
#'   train <- data[1:150, ]
#'   test <- data[151:nrow(data), ]
#'   fit <- brms::brm(
#'     Reaction ~ Days, data = train,
#'     refresh = 0, chains = 2, iter = 1000
#'   )
#'   test_pred_measure(
#'     y = test$Reaction,
#'     ypred = brms::posterior_predict(fit, newdata = test),
#'     mupred = brms::posterior_epred(fit, newdata = test),
#'     ylp = brms::log_lik(fit),
#'     ylp_test = brms::log_lik(fit, newdata = test),
#'     measure = c("rmse", "r2")
#'   )
#' }
#' }
#'
#' @seealso [insample_pred_measure()], [loo_pred_measure()],
#'   [kfold_pred_measure()], [pred_measure()], [supported_measures_list],
#'   [pred-measure workflow article](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html)
#'
#' @export
test_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_test = NULL,
  measure = NULL,
  group_ids = NULL,
  control = list()
) {
  do_pred_measure(
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_test = ylp_test,
    measure = measure, 
    predperf = NULL,
    loo = NULL, 
    kfold = NULL,
    group_ids = group_ids,
    psis_object = NULL,
    save_psis = FALSE,
    source = "test",
    control = control
  )
}

#' Add predictive performance measures to an existing result
#'
#' @description
#' Extend a `"pred_measure"` object with additional measures **without
#' recomputing** what is already stored. Use this for interactive exploration
#' or when you first compute base density summaries and later add distributional
#' or point-prediction metrics.
#'
#' Pass the existing object as `predperf` and supply any inputs newly required
#' by the requested measures (see the input table in
#' [insample_pred_measure()]). The evaluation mode (`"insample"`, `"loo"`,
#' `"kfold"`, or `"test"`) is taken from `predperf`; LOO paths reuse stored
#' PSIS weights when available.
#'
#' @inheritParams pred_measure_params
#'
#' @return
#' An updated object of the same class as `predperf`, with new rows in
#' `estimates` and columns in `pointwise` for each requested measure. Base
#' summaries (`elpd` and LOO/k-fold complexity terms such as `p_loo`) are not
#' recomputed.
#'
#' @details
#' **Typical workflow:**
#'
#' \preformatted{
#' result <- loo_pred_measure(loo = loo_fit, save_psis = TRUE)
#' pred_measure(
#'   y = y,
#'   mupred = mupred,
#'   predperf = result,
#'   measure = c("rmse", "r2")
#' )
#' }
#'
#' When extending a LOO result, ensure the initial call used `save_psis = TRUE`
#' (or that `predperf` already contains a `psis_object`) so LOO weights are
#' available for additional measures.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit <- brms::brm(
#'     Reaction ~ Days, data = lme4::sleepstudy,
#'     refresh = 0, chains = 2, iter = 1000
#'   )
#'   result <- insample_pred_measure(
#'     ylp = brms::log_lik(fit),
#'     y = fit$data$Reaction,
#'     ypred = brms::posterior_predict(fit),
#'     measure = "rmse"
#'   )
#'   pred_measure(
#'     y = fit$data$Reaction,
#'     mupred = brms::posterior_epred(fit),
#'     predperf = result,
#'     measure = "r2"
#'   )
#' }
#' }
#'
#' @seealso [insample_pred_measure()], [loo_pred_measure()],
#'   [kfold_pred_measure()], [test_pred_measure()], [supported_measures_list],
#'   [pred-measure workflow article](https://mc-stan.org/loo/articles/articles-online-only/pred-measure-workflow.html)
#'
#' @export
pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  measure = NULL,
  predperf,
  group_ids = NULL,
  psis_object = NULL,
  save_psis = FALSE,
  control = list()
) {
  do_pred_measure(
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_test = NULL,
    measure = measure, 
    predperf = predperf,
    loo = NULL,
    kfold = NULL,
    group_ids = group_ids,
    psis_object = psis_object,
    save_psis = save_psis,
    source = attr(predperf, "source"),
    control = control
  )
}