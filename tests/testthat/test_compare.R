set.seed(123)

LLarr <- example_loglik_array()
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- suppressWarnings(waic(LLarr))
w2 <- suppressWarnings(waic(LLarr2))

test_that("loo_compare throws appropriate errors", {
  w3 <- suppressWarnings(waic(LLarr[,, -1]))
  w4 <- suppressWarnings(waic(LLarr[,, -(1:2)]))

  expect_error(loo_compare(2, 3), "must be a list if not a 'loo' or 'pred_measure' object")
  expect_error(
    loo_compare(w1, w2, x = list(w1, w2)),
    "If 'x' is a list then '...' should not be specified"
  )
  expect_error(loo_compare(w1, list(1, 2, 3)), "class 'loo'")
  expect_error(loo_compare(w1), "requires at least two models")
  expect_error(loo_compare(x = list(w1)), "requires at least two models")
  expect_error(
    loo_compare(w1, w3),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'model1' (32), 'model2' (31)",
    fixed = TRUE
  )
  expect_error(
    loo_compare(w1, w2, w3),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'model1' (32), 'model2' (32), 'model3' (31)",
    fixed = TRUE
  )
  expect_error(
    loo_compare(x = list("Model A" = w1, "Model B" = w2, "Model C" = w3)),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'Model A' (32), 'Model B' (32), 'Model C' (31)",
    fixed = TRUE
  )
})

test_that("loo_compare dispatches loo_pred_measure inputs", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mse")
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mse")
  )

  comp <- suppressMessages(loo_compare(pm1, pm2))
  expect_s3_class(comp, "compare.loo")
  expect_null(attr(comp, "rank_by"))
  expect_true(all(c("elpd_diff", "se_diff", "p_worse", "diag_diff") %in% colnames(comp)))
  expect_true(all(c("r2_diff", "r2_se_diff", "mse_diff", "mse_se_diff") %in% colnames(comp)))
  expect_false(anyNA(comp$r2_se_diff))
  expect_false(anyNA(comp$mse_se_diff))
  expect_false("r2_loo_diff" %in% colnames(comp))
  expect_false("mse_p_worse" %in% colnames(comp))

  expect_error(
    loo_compare(w1, pm1),
    "Cannot mix 'loo_pred_measure' objects with other 'loo' objects",
    fixed = TRUE
  )
  expect_error(
    loo_compare(
      insample_pred_measure(ylp = res$ylp_m1),
      insample_pred_measure(ylp = res$ylp_m1)
    ),
    "requires 'loo_pred_measure' objects",
    fixed = TRUE
  )
  expect_error(loo_compare(pm1), "requires at least two models", fixed = TRUE)
  expect_null(attr(loo_compare(w1, w2), "rank_by"))
})

test_that("loo_compare warns when predictive measures differ across models", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mse")
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = c("r2", "mae")
  )

  expect_warning(
    comp <- suppressMessages(loo_compare(list(m1 = pm1, m2 = pm2))),
    "Omitted measures: mae \\(m2\\), mse \\(m1\\)"
  )
  expect_equal(attr(comp, "compare_measures"), c("elpd", "r2"))
  expect_false("mse_diff" %in% colnames(comp))
  expect_false("mae_diff" %in% colnames(comp))
})

test_that("loo_compare works with three loo_pred_measure models", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mae")
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = c("r2", "mae")
  )
  pm3 <- loo_pred_measure(
    loo = res$loo_p_m3,
    y = res$y,
    mupred = res$mupred_m3,
    ylp = res$ylp_m3,
    measure = c("r2", "mae")
  )

  comp <- loo_compare(
    list("A" = pm1, "B" = pm2, "C" = pm3),
    rank_by = "mae"
  )
  expect_snapshot(print(comp))
  expect_equal(nrow(comp), 3L)
  expect_equal(comp$model, c("C", "B", "A"))
  expect_equal(attr(comp, "rank_by"), "mae")
  expect_equal(attr(comp, "compare_measures"), c("elpd", "r2", "mae"))
  expect_equal(comp$mae_diff[1L], 0)
  expect_true(all(comp$mae_diff[-1L] < 0))
  expect_true(all(comp$elpd_diff[-1L] <= 0))
  expect_equal(attr(comp, "sign_converted_measures"), c("mae"))
})

test_that("loo_compare informs when measure signs are converted", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mse")
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = c("r2", "mse")
  )

  expect_snapshot(comp <- loo_compare(pm1, pm2))
  expect_equal(attr(comp, "sign_converted_measures"), "mse")

  pm_elpd <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1
  )
  expect_no_message(loo_compare(pm_elpd, pm_elpd))
})

test_that("loo_compare rank_by changes order for loo_pred_measure", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mae")
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = c("r2", "mae")
  )

  comp_elpd <- loo_compare(pm1, pm2, rank_by = "elpd")
  comp_mse <- loo_compare(pm1, pm2, rank_by = "mae")
  expect_equal(attr(comp_elpd, "rank_by"), "elpd")
  expect_equal(attr(comp_mse, "rank_by"), "mae")
  expect_equal(comp_elpd$elpd_diff[1L], 0)
  expect_equal(comp_mse$mae_diff[1L], 0)
})

test_that("print.compare.loo works for loo_pred_measure comparisons", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mae")
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = c("r2", "mae")
  )
  pm3 <- loo_pred_measure(
    loo = res$loo_p_m3,
    y = res$y,
    mupred = res$mupred_m3,
    ylp = res$ylp_m3,
    measure = c("r2", "mae")
  )

  comp <- suppressMessages(loo_compare(list(m1 = pm1, m2 = pm2, m3 = pm3)))
  expect_snapshot(print(comp))
  expect_snapshot(print(comp, measures = "all", digits = 2))
  expect_snapshot(print(comp, measures = c("r2", "mae")))

  comp_mae <- suppressMessages(loo_compare(list(m1 = pm1, m2 = pm2), rank_by = "mae"))
  expect_snapshot(print(comp_mae))

  expect_error(
    print(comp, measures = "foo"),
    "Unknown measure\\(s\\) in `measures`"
  )
})

test_that("loo_compare measure helpers work as expected", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = c("r2", "mse")
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = c("r2", "mse")
  )
  loos <- list(pm1, pm2)
  cols <- loo:::.compare_pointwise_cols(loos)

  expect_equal(cols, c("elpd_loo", "r2_loo", "mse_loo"))
  expect_equal(loo:::.compare_measures(loos), c("elpd", "r2", "mse"))
  expect_equal(loo:::.pointwise_col("mse", cols), "mse_loo")
  expect_equal(loo:::.pointwise_col("elpd", cols), "elpd_loo")
  expect_equal(loo:::.display_name("rmse_loo"), "rmse")
  expect_equal(loo:::.resolve_rank_measure(loos, NULL)$bare, "elpd")
  expect_equal(loo:::.resolve_rank_measure(loos, "mse")$internal, "mse_loo")
  expect_true(loo:::.is_elpd_measure("elpd_loo"))
  expect_false(loo:::.is_elpd_measure("mse_loo"))
  expect_equal(attr(pm1, "measure_higher_is_better")$mse, NULL)
  expect_equal(attr(pm1, "measure_higher_is_better")$r2, NULL)
  expect_equal(attr(pm1, "measure_higher_is_better")$elpd, NULL)
  expect_equal(attr(pm1, "measure_compare_meta")$elpd$diff_method, "sum")
  expect_equal(attr(pm1, "measure_compare_meta")$mse$loss, TRUE)
  expect_equal(attr(pm1, "measure_compare_meta")$mse$diff_method, "mean")
  expect_equal(attr(pm1, "measure_compare_meta")$r2$diff_method, "pairwise")
  expect_equal(attr(pm1, "measure_compare_meta")$r2$se_diff_fun, "r2")
  expect_equal(
    attr(pm1, "measure_compare_meta")$r2$extra$mse_y_i,
    (res$y - mean(res$y))^2
  )
  # only measures that need it carry `extra`
  expect_null(attr(pm1, "measure_compare_meta")$mse$extra)
  expect_true(loo:::.measure_lower_is_better("mse_loo", loos))
  expect_false(loo:::.measure_lower_is_better("r2_loo", loos))
  expect_true(loo:::.measure_lower_is_better("mse_loo"))
  expect_false(loo:::.measure_lower_is_better("r2_loo"))
  expect_equal(
    loo:::.compare_sign_converted_measures(c("elpd_loo", "mse_loo", "r2_loo"), loos),
    c("mse")
  )

  pair_stats_elpd <- loo:::.pair_measure_stats(
    pm2, pm1, "elpd_loo", "sum", loos = loos
  )
  expect_equal(unname(pair_stats_elpd["se"]), loo:::se_elpd_diff(
    pm2$pointwise[, "elpd_loo"] - pm1$pointwise[, "elpd_loo"]
  ))
  expect_equal(
    unname(loo:::.pair_measure_stats(pm1, pm1, "elpd_loo", "sum", loos = loos)["diff"]),
    0
  )

  pair_mse <- loo:::.pair_measure_stats(
    pm2, pm1, "mse_loo", "mean", loos = loos
  )
  expect_equal(
    unname(pair_mse["diff"]),
    pm1$estimates["mse_loo", "Estimate"] - pm2$estimates["mse_loo", "Estimate"]
  )
  pair_r2 <- loo:::.pair_measure_stats(
    pm2, pm1, "r2_loo", "pairwise", loos = loos
  )
  expect_equal(
    unname(pair_r2["diff"]),
    pm2$estimates["r2_loo", "Estimate"] - pm1$estimates["r2_loo", "Estimate"]
  )
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "mse_loo"), "mean")
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "r2_loo"), "pairwise")
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "elpd_loo"), "sum")
  expect_equal(
    unname(pair_mse["se"]),
    stats::sd(
      pm2$pointwise[, "mse_loo"] - pm1$pointwise[, "mse_loo"]
    ) / sqrt(nrow(pm1$pointwise))
  )
})

test_that("rmse differences use the delta-method standard error", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = "rmse"
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = "rmse"
  )
  loos <- list(pm1, pm2)

  expect_equal(loo:::.measure_pointwise_diff_method(loos, "rmse_loo"), "pairwise")

  pair <- loo:::.pair_measure_stats(pm2, pm1, "rmse_loo", loos = loos)

  # rmse is a loss, so the reported difference is on the utility scale
  expect_equal(
    unname(pair["diff"]),
    pm1$estimates["rmse_loo", "Estimate"] - pm2$estimates["rmse_loo", "Estimate"]
  )

  # first-order bivariate Taylor approximation propagated from the MSE scale,
  # using the covariance between the two models' pointwise squared errors
  sqe1 <- pm1$pointwise[, "rmse_loo"]
  sqe2 <- pm2$pointwise[, "rmse_loo"]
  n <- length(sqe1)
  mse1 <- mean(sqe1)
  mse2 <- mean(sqe2)
  cov_mse <- sum((sqe2 - mse2) * (sqe1 - mse1)) / (n * (n - 1))
  expected_se <- 0.5 * sqrt(
    (var(sqe2) / n) / mse2 +
      (var(sqe1) / n) / mse1 -
      2 * cov_mse / sqrt(mse2 * mse1)
  )
  expect_equal(unname(pair["se"]), expected_se)

  # the standard error is a proper paired quantity, not a sum of the two
  # per-model standard errors
  expect_lt(
    unname(pair["se"]),
    pm1$estimates["rmse_loo", "SE"] + pm2$estimates["rmse_loo", "SE"]
  )

  # a model compared against itself has zero difference and zero uncertainty
  self <- loo:::.pair_measure_stats(pm1, pm1, "rmse_loo", loos = loos)
  expect_equal(unname(self["diff"]), 0)
  expect_equal(unname(self["se"]), 0)

  comp <- suppressMessages(loo_compare(pm1, pm2))
  expect_false(anyNA(comp$rmse_se_diff))
  expect_equal(comp$rmse_se_diff[1], 0)
})

test_that("r2 differences use the delta-method standard error", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  make <- function(loo, mupred, ylp) {
    loo_pred_measure(
      loo = loo, y = res$y, mupred = mupred, ylp = ylp, measure = "r2"
    )
  }
  pm1 <- make(res$loo_p_m1, res$mupred_m1, res$ylp_m1)
  pm2 <- make(res$loo_p_m2, res$mupred_m2, res$ylp_m2)
  loos <- list(pm1, pm2)

  expect_equal(loo:::.measure_pointwise_diff_method(loos, "r2_loo"), "pairwise")

  pair <- loo:::.pair_measure_stats(pm2, pm1, "r2_loo", loos = loos)

  # r2 is already a utility, so the difference is reported as stored
  expect_equal(
    unname(pair["diff"]),
    pm2$estimates["r2_loo", "Estimate"] - pm1$estimates["r2_loo", "Estimate"]
  )

  # first-order trivariate Taylor approximation, written out term by term as in
  # the derivation rather than in the collapsed single-variance form the
  # implementation uses
  sqe1 <- pm1$pointwise[, "r2_loo"]
  sqe2 <- pm2$pointwise[, "r2_loo"]
  d <- sqe2 - sqe1
  n <- length(d)
  mse_diff <- mean(d)
  msey_i <- (res$y - mean(res$y))^2
  mse_y <- mean(msey_i)
  t1 <- var(d) / n
  t2 <- -2 * (mse_diff / mse_y) *
    (sum((d - mse_diff) * (msey_i - mse_y)) / (n * (n - 1)))
  t3 <- (mse_diff^2 / mse_y^2) * (var(msey_i) / n)
  expect_equal(unname(pair["se"]), sqrt(t1 + t2 + t3) / mse_y)

  # the difference is also the negative MSE difference over the baseline
  expect_equal(unname(pair["diff"]), -mse_diff / mse_y)

  # the uncertainty in a difference does not depend on which model is the
  # reference, even though the difference itself changes sign
  swapped <- loo:::.pair_measure_stats(pm1, pm2, "r2_loo", loos = loos)
  expect_equal(unname(swapped["se"]), unname(pair["se"]))
  expect_equal(unname(swapped["diff"]), -unname(pair["diff"]))

  # a model compared against itself has zero difference and zero uncertainty
  self <- loo:::.pair_measure_stats(pm1, pm1, "r2_loo", loos = loos)
  expect_equal(unname(self["diff"]), 0)
  expect_equal(unname(self["se"]), 0)

  comp <- suppressMessages(loo_compare(pm1, pm2))
  expect_false(anyNA(comp$r2_se_diff))
  expect_equal(comp$r2_se_diff[1], 0)

  # the single-model standard error is the same expansion evaluated at one
  # model's squared errors
  t1 <- var(sqe1) / n
  t2 <- -2 * (mean(sqe1) / mse_y) *
    (sum((sqe1 - mean(sqe1)) * (msey_i - mse_y)) / (n * (n - 1)))
  t3 <- (mean(sqe1)^2 / mse_y^2) * (var(msey_i) / n)
  expect_equal(
    unname(pm1$estimates["r2_loo", "SE"]),
    sqrt(t1 + t2 + t3) / mse_y
  )
})

test_that("r2 se_diff is invariant to `higher_is_better`", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  make <- function(loo, mupred, ylp, higher_is_better) {
    loo_pred_measure(
      loo = loo,
      y = res$y,
      mupred = mupred,
      ylp = ylp,
      measure = "r2",
      control = list(r2 = list(higher_is_better = higher_is_better))
    )
  }

  utility <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, TRUE),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, TRUE)
  )
  loss <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, FALSE),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, FALSE)
  )

  # r2 is naturally a utility, so `higher_is_better = FALSE` is the flipped one
  expect_equal(loo:::.measure_natural_sign("r2_loo", utility), 1)
  expect_equal(loo:::.measure_natural_sign("r2_loo", loss), -1)

  # the baseline is a property of `y`, so it is never sign-flipped
  expect_equal(
    attr(loss[[1L]], "measure_compare_meta")$r2$extra$mse_y_i,
    attr(utility[[1L]], "measure_compare_meta")$r2$extra$mse_y_i
  )

  pair_utility <- loo:::.pair_measure_stats(
    utility[[2L]], utility[[1L]], "r2_loo", loos = utility
  )
  pair_loss <- loo:::.pair_measure_stats(
    loss[[2L]], loss[[1L]], "r2_loo", loos = loss
  )
  expect_equal(unname(pair_loss["se"]), unname(pair_utility["se"]))
  expect_equal(unname(pair_loss["diff"]), unname(pair_utility["diff"]))
})

test_that("r2 reports the difference without an se when the baseline is gone", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  make <- function(loo, mupred, ylp) {
    loo_pred_measure(
      loo = loo, y = res$y, mupred = mupred, ylp = ylp, measure = "r2"
    )
  }
  pm1 <- make(res$loo_p_m1, res$mupred_m1, res$ylp_m1)
  pm2 <- make(res$loo_p_m2, res$mupred_m2, res$ylp_m2)

  # objects computed before the baseline was stored
  drop_baseline <- function(x) {
    meta <- attr(x, "measure_compare_meta")
    meta$r2$extra <- NULL
    attr(x, "measure_compare_meta") <- meta
    x
  }
  stale1 <- drop_baseline(pm1)
  stale2 <- drop_baseline(pm2)

  pair <- loo:::.pair_measure_stats(
    stale2, stale1, "r2_loo", loos = list(stale1, stale2)
  )
  expect_equal(
    unname(pair["diff"]),
    pm2$estimates["r2_loo", "Estimate"] - pm1$estimates["r2_loo", "Estimate"]
  )
  expect_true(is.na(pair["se"]))

  comp <- suppressMessages(loo_compare(stale1, stale2))
  expect_false(anyNA(comp$r2_diff))
  expect_true(all(is.na(comp$r2_se_diff)))

  # one stale model does not cost the others their standard error: the
  # baseline is shared, so the other model's copy is used, and the metadata
  # check ignores `extra` rather than reporting it as a disagreement
  mixed <- suppressMessages(loo_compare(stale1, pm2))
  expect_false(anyNA(mixed$r2_se_diff))
})

test_that("rmse se_diff is invariant to `higher_is_better`", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  make <- function(loo, mupred, ylp, higher_is_better) {
    loo_pred_measure(
      loo = loo,
      y = res$y,
      mupred = mupred,
      ylp = ylp,
      measure = "rmse",
      control = list(rmse = list(higher_is_better = higher_is_better))
    )
  }

  loss <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, FALSE),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, FALSE)
  )
  utility <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, TRUE),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, TRUE)
  )

  # stored on opposite scales, so the natural scale must be restored before the
  # square roots and ratios of the delta method are applied
  expect_equal(loo:::.measure_natural_sign("rmse_loo", loss), 1)
  expect_equal(loo:::.measure_natural_sign("rmse_loo", utility), -1)

  pair_loss <- loo:::.pair_measure_stats(
    loss[[2L]], loss[[1L]], "rmse_loo", loos = loss
  )
  pair_utility <- loo:::.pair_measure_stats(
    utility[[2L]], utility[[1L]], "rmse_loo", loos = utility
  )
  expect_equal(unname(pair_loss["se"]), unname(pair_utility["se"]))
  expect_equal(unname(pair_loss["diff"]), unname(pair_utility["diff"]))
})

# two balanced-accuracy measures over the same three-class outcome: the second
# model has probability mass shifted towards the first (and largest) class, so
# the two disagree on a subset of observations and the class strata are
# unbalanced, which is where balanced accuracy differs from plain accuracy
.make_bacc_pms <- function(bias = 0.6, higher_is_better = NULL) {
  res <- readRDS("data-for-tests/test_data_penguins.Rds")
  y <- as.integer(res$y)
  set.seed(4321)
  ylp <- matrix(
    rnorm(nrow(res$mupred) * ncol(res$mupred)),
    nrow = nrow(res$mupred)
  )
  biased <- res$mupred
  biased[, , 1L] <- biased[, , 1L] + bias
  biased <- sweep(biased, c(1, 2), apply(biased, c(1, 2), sum), "/")

  make <- function(mupred) {
    suppressWarnings(loo_pred_measure(
      ylp = ylp,
      y = y,
      mupred = mupred,
      measure = "bacc",
      control = list(bacc = list(higher_is_better = higher_is_better))
    ))
  }
  list(pm1 = make(res$mupred), pm2 = make(biased), y = y)
}

test_that("bacc differences use the stratified paired standard error", {
  fx <- .make_bacc_pms()
  pm1 <- fx$pm1
  pm2 <- fx$pm2
  loos <- list(pm1, pm2)

  expect_equal(loo:::.measure_pointwise_diff_method(loos, "bacc_loo"), "pairwise")
  expect_equal(
    attr(pm1, "measure_compare_meta")$bacc$se_diff_fun,
    "bacc"
  )

  pair <- loo:::.pair_measure_stats(pm2, pm1, "bacc_loo", loos = loos)

  # bacc is already a utility, so the difference is reported as stored
  expect_equal(
    unname(pair["diff"]),
    pm2$estimates["bacc_loo", "Estimate"] - pm1$estimates["bacc_loo", "Estimate"]
  )

  # recover the 0/1 accuracies and check the standard error against the
  # McNemar discordant-count form of the paired difference of proportions,
  # written out per stratum rather than in the pointwise-variance form the
  # implementation uses
  class_id <- attr(pm1, "measure_compare_meta")$bacc$extra$class_id
  n_c <- tabulate(class_id)
  K <- length(n_c)
  acc1 <- round(pm1$pointwise[, "bacc_loo"] * K * n_c[class_id])
  acc2 <- round(pm2$pointwise[, "bacc_loo"] * K * n_c[class_id])
  expect_true(all(acc1 %in% c(0, 1)) && all(acc2 %in% c(0, 1)))

  v <- 0
  for (k in seq_len(K)) {
    in_k <- class_id == k
    b <- sum(acc2[in_k] == 1 & acc1[in_k] == 0)
    cc <- sum(acc2[in_k] == 0 & acc1[in_k] == 1)
    nk <- n_c[k]
    # the Wald paired-proportion variance, scaled by nk / (nk - 1) to match the
    # sample variance the implementation takes
    v <- v + ((b + cc) / nk^2 - (b - cc)^2 / nk^3) * (nk / (nk - 1))
  }
  expect_equal(unname(pair["se"]), sqrt(v) / K)

  # the strata carry information the pointwise vector alone does not: pooling
  # them would give a different answer, so this is not the `"mean"` path
  d <- (pm2$pointwise[, "bacc_loo"] - pm1$pointwise[, "bacc_loo"])
  expect_false(isTRUE(all.equal(
    unname(pair["se"]),
    sd(d) / sqrt(length(d))
  )))

  # the uncertainty in a difference does not depend on which model is the
  # reference, even though the difference itself changes sign
  swapped <- loo:::.pair_measure_stats(pm1, pm2, "bacc_loo", loos = loos)
  expect_equal(unname(swapped["se"]), unname(pair["se"]))
  expect_equal(unname(swapped["diff"]), -unname(pair["diff"]))

  # a model compared against itself has zero difference and zero uncertainty
  self <- loo:::.pair_measure_stats(pm1, pm1, "bacc_loo", loos = loos)
  expect_equal(unname(self["diff"]), 0)
  expect_equal(unname(self["se"]), 0)

  comp <- suppressMessages(loo_compare(pm1, pm2))
  expect_false(anyNA(comp$bacc_se_diff))
  expect_equal(comp$bacc_se_diff[1], 0)
})

test_that("bacc se_diff is invariant to `higher_is_better`", {
  utility <- .make_bacc_pms(higher_is_better = TRUE)
  loss <- .make_bacc_pms(higher_is_better = FALSE)

  # bacc is naturally a utility, so `higher_is_better = FALSE` is the flipped one
  expect_equal(
    loo:::.measure_natural_sign("bacc_loo", list(utility$pm1, utility$pm2)), 1
  )
  expect_equal(
    loo:::.measure_natural_sign("bacc_loo", list(loss$pm1, loss$pm2)), -1
  )

  # the strata are a property of `y`, so they are never sign-flipped
  expect_equal(
    attr(loss$pm1, "measure_compare_meta")$bacc$extra$class_id,
    attr(utility$pm1, "measure_compare_meta")$bacc$extra$class_id
  )

  pair_utility <- loo:::.pair_measure_stats(
    utility$pm2, utility$pm1, "bacc_loo", loos = list(utility$pm1, utility$pm2)
  )
  pair_loss <- loo:::.pair_measure_stats(
    loss$pm2, loss$pm1, "bacc_loo", loos = list(loss$pm1, loss$pm2)
  )
  expect_equal(unname(pair_loss["se"]), unname(pair_utility["se"]))
  expect_equal(unname(pair_loss["diff"]), unname(pair_utility["diff"]))
})

test_that("bacc reports the difference without an se when the strata are gone", {
  fx <- .make_bacc_pms()

  # objects computed before the class strata were stored
  drop_strata <- function(x) {
    meta <- attr(x, "measure_compare_meta")
    meta$bacc$extra <- NULL
    attr(x, "measure_compare_meta") <- meta
    x
  }
  stale1 <- drop_strata(fx$pm1)
  stale2 <- drop_strata(fx$pm2)

  pair <- loo:::.pair_measure_stats(
    stale2, stale1, "bacc_loo", loos = list(stale1, stale2)
  )
  expect_equal(
    unname(pair["diff"]),
    fx$pm2$estimates["bacc_loo", "Estimate"] -
      fx$pm1$estimates["bacc_loo", "Estimate"]
  )
  expect_true(is.na(pair["se"]))

  # the strata are shared, so a stale model paired with a current one still
  # gets a standard error, from whichever copy survives. `pm1` is the better
  # model and so heads the table; staleness in `pm2` costs nothing at all
  mixed <- suppressMessages(loo_compare(fx$pm1, stale2))
  expect_false(anyNA(mixed$bacc_se_diff))

  # but a stale model at the head of the table has no second copy to fall back
  # on for its own row, which is a comparison against itself
  mixed_stale_first <- suppressMessages(loo_compare(stale1, fx$pm2))
  expect_true(is.na(mixed_stale_first$bacc_se_diff[1L]))
  expect_false(is.na(mixed_stale_first$bacc_se_diff[2L]))
})

test_that("custom measures can supply their own se_diff_fun", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")

  # a custom rmse, whose estimate is neither a sum nor a mean of `pointwise`
  my_rmse <- function(y, mupred) {
    sqe_i <- (y - colMeans(mupred))^2
    list(
      estimate = sqrt(mean(sqe_i)),
      se = sqrt(var(sqe_i) / length(sqe_i)) / (2 * sqrt(mean(sqe_i))),
      pointwise = sqe_i
    )
  }
  attr(my_rmse, "measure_name") <- "my_rmse"

  make <- function(loo, mupred, ylp, fun) {
    loo_pred_measure(
      loo = loo, y = res$y, mupred = mupred, ylp = ylp, measure = fun
    )
  }

  # without an se_diff_fun the difference has no standard error
  plain <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_rmse),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_rmse)
  )
  expect_equal(
    loo:::.measure_pointwise_diff_method(plain, "my_rmse_loo"),
    "estimates_only"
  )
  expect_true(is.na(
    loo:::.pair_measure_stats(plain[[2L]], plain[[1L]], "my_rmse_loo", loos = plain)["se"]
  ))

  attr(my_rmse, "se_diff_fun") <- loo:::.se_diff_rmse
  with_se <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_rmse),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_rmse)
  )
  expect_equal(
    attr(with_se[[1L]], "measure_compare_meta")$my_rmse$diff_method,
    "pairwise"
  )
  expect_equal(
    loo:::.measure_pointwise_diff_method(with_se, "my_rmse_loo"),
    "pairwise"
  )

  pair <- loo:::.pair_measure_stats(
    with_se[[2L]], with_se[[1L]], "my_rmse_loo", loos = with_se
  )
  expect_false(is.na(pair["se"]))
  expect_gt(unname(pair["se"]), 0)

  # a custom measure can carry its own auxiliary data through to its
  # se_diff_fun, and each model receives its own copy
  my_scaled <- function(y, mupred) {
    sqe_i <- (y - colMeans(mupred))^2
    list(
      estimate = sqrt(mean(sqe_i)),
      se = sqrt(var(sqe_i) / length(sqe_i)) / (2 * sqrt(mean(sqe_i))),
      pointwise = sqe_i,
      extra = list(scale = stats::sd(y), n_used = length(y))
    )
  }
  attr(my_scaled, "measure_name") <- "my_scaled"
  attr(my_scaled, "se_diff_fun") <- function(ref, cmp) {
    stopifnot(
      identical(ref$extra$n_used, length(ref$pointwise)),
      identical(cmp$extra$scale, ref$extra$scale)
    )
    loo:::.se_diff_rmse(ref, cmp) / ref$extra$scale
  }
  scaled <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_scaled),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_scaled)
  )
  expect_equal(
    attr(scaled[[1L]], "measure_compare_meta")$my_scaled$extra,
    list(scale = stats::sd(res$y), n_used = length(res$y))
  )
  pair_scaled <- loo:::.pair_measure_stats(
    scaled[[2L]], scaled[[1L]], "my_scaled_loo", loos = scaled
  )
  expect_equal(
    unname(pair_scaled["se"]),
    unname(pair["se"]) / stats::sd(res$y)
  )

  # `extra` that is not a list is rejected at compute time
  my_bad_extra <- function(y, mupred) {
    sqe_i <- (y - colMeans(mupred))^2
    list(
      estimate = mean(sqe_i),
      se = sqrt(var(sqe_i) / length(sqe_i)),
      pointwise = sqe_i,
      extra = 1
    )
  }
  attr(my_bad_extra, "measure_name") <- "my_bad_extra"
  expect_error(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_bad_extra),
    "must be a list"
  )

  # a custom measure declaring an se_diff_fun that returns nonsense is caught
  attr(my_rmse, "se_diff_fun") <- function(ref, cmp) c(1, 2)
  bad <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_rmse),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_rmse)
  )
  expect_error(
    loo:::.pair_measure_stats(bad[[2L]], bad[[1L]], "my_rmse_loo", loos = bad),
    "must return a numeric scalar"
  )
})

test_that("loo_compare errors on inconsistent measure metadata", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = "mse",
    control = list(mse = list(higher_is_better = NULL))
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = "mse",
    control = list(mse = list(higher_is_better = TRUE))
  )

  expect_error(
    suppressMessages(loo_compare(pm1, pm2)),
    "disagree on comparison metadata for measure 'mse'"
  )
})

test_that("loo_compare errors when compare metadata is missing on some models", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  pm1 <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1,
    measure = "mse"
  )
  pm2 <- loo_pred_measure(
    loo = res$loo_p_m2,
    y = res$y,
    mupred = res$mupred_m2,
    ylp = res$ylp_m2,
    measure = "mse"
  )
  compare_meta <- attr(pm2, "measure_compare_meta")
  compare_meta$mse <- NULL
  attr(pm2, "measure_compare_meta") <- compare_meta

  expect_error(
    suppressMessages(loo_compare(pm1, pm2)),
    "Not all models provide comparison metadata for measure 'mse'"
  )
})

test_that("loo_compare warns when rank_by is ignored for classic loo objects", {
  expect_warning(
    loo_compare(w1, w2, rank_by = "mse"),
    "`rank_by` is only used for `loo_pred_measure` comparisons"
  )
})

.make_compare_pm <- function(res, model = 1L, measure, extra_args = list()) {
  suffix <- model
  args <- c(
    list(
      loo = res[[paste0("loo_p_m", suffix)]],
      y = res$y,
      mupred = res[[paste0("mupred_m", suffix)]],
      ylp = res[[paste0("ylp_m", suffix)]],
      measure = measure
    ),
    extra_args
  )
  do.call(loo_pred_measure, args)
}

.make_compare_pm_synthetic <- function(measure) {
  if (measure == "brier") {
    res_binary <- readRDS("data-for-tests/test_data_binary.Rds")
    ylp <- matrix(
      rnorm(nrow(res_binary$ypred) * ncol(res_binary$ypred)),
      nrow = nrow(res_binary$ypred)
    )
    return(loo_pred_measure(
      ylp = ylp,
      y = res_binary$y,
      ypred = res_binary$ypred,
      measure = measure
    ))
  }
  if (measure %in% c("acc", "bacc")) {
    res_cat <- readRDS("data-for-tests/test_data_penguins.Rds")
    ylp <- matrix(
      rnorm(nrow(res_cat$mupred) * ncol(res_cat$mupred)),
      nrow = nrow(res_cat$mupred)
    )
    return(loo_pred_measure(
      ylp = ylp,
      y = as.integer(res_cat$y),
      mupred = res_cat$mupred,
      measure = measure
    ))
  }
  stop("Unsupported synthetic measure: ", measure)
}

test_that("loo_compare works for all built-in measures", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  res_roaches <- readRDS("data-for-tests/test_data_roaches.Rds")
  roaches_measures <- c(
    "ic", "mlpd", "mae", "r2", "rmse", "mse"
  )
  for (measure in roaches_measures) {
    pm1 <- .make_compare_pm(res, 1L, measure)
    pm2 <- .make_compare_pm(res, 2L, measure)
    comp <- suppressMessages(loo_compare(pm1, pm2))
    expect_true(paste0(measure, "_diff") %in% colnames(comp), info = measure)
    expect_equal(attr(comp, "compare_measures"), c("elpd", measure), info = measure)
  }

  for (measure in c("rps", "srps")) {
    pm1 <- loo_pred_measure(
      loo = res$loo_p_m1,
      y = res$y,
      ypred = res_roaches$ypred,
      ylp = res$ylp_m1,
      measure = measure
    )
    pm2 <- loo_pred_measure(
      loo = res$loo_p_m2,
      y = res$y,
      ypred = res_roaches$ypred,
      ylp = res$ylp_m2,
      measure = measure
    )
    comp <- suppressMessages(loo_compare(pm1, pm2))
    expect_true(paste0(measure, "_diff") %in% colnames(comp), info = measure)
    expect_equal(attr(comp, "compare_measures"), c("elpd", measure), info = measure)
  }

  for (measure in c("brier", "acc", "bacc")) {
    pm1 <- .make_compare_pm_synthetic(measure)
    pm2 <- .make_compare_pm_synthetic(measure)
    comp <- suppressMessages(loo_compare(pm1, pm2))
    expect_true(paste0(measure, "_diff") %in% colnames(comp), info = measure)
    expect_equal(attr(comp, "compare_measures"), c("elpd", measure), info = measure)
  }
})

.make_many_compare_pms <- function(res, n, noise_scale = 0.01) {
  lapply(seq_len(n), function(i) {
    loo_pred_measure(
      loo = res$loo_p_m1,
      y = res$y,
      mupred = res$mupred_m1 + rnorm(length(res$y), 0, noise_scale * i),
      ylp = res$ylp_m1,
      measure = "mae"
    )
  })
}

test_that("loo_compare warns for many loo_pred_measure models", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  set.seed(123)
  pm_list <- .make_many_compare_pms(res, 25L)
  expect_warning(
    suppressMessages(loo_compare(pm_list)),
    "Difference in performance potentially due to chance. See McLatchie and Vehtari (2023) for details.",
    fixed = TRUE
  )

  pm_list_short <- .make_many_compare_pms(res, 4L)
  expect_no_warning(suppressMessages(loo_compare(pm_list_short)))
})

test_that("loo_compare throws appropriate warnings", {
  w3 <- w1
  w4 <- w2
  class(w3) <- class(w4) <- c("kfold", "loo")
  attr(w3, "K") <- 2
  attr(w4, "K") <- 3
  expect_warning(
    loo_compare(w3, w4),
    "Not all kfold objects have the same K value"
  )

  class(w4) <- c("psis_loo", "loo")
  attr(w4, "K") <- NULL
  expect_warning(loo_compare(w3, w4), "Comparing LOO-CV to K-fold-CV")

  w3 <- w1
  w4 <- w2
  attr(w3, "yhash") <- "a"
  attr(w4, "yhash") <- "b"
  expect_warning(loo_compare(w3, w4), "Not all models have the same y variable")

  set.seed(123)
  w_list <- lapply(1:25, function(x) {
    suppressWarnings(waic(LLarr + rnorm(1, 0, 0.1)))
  })
  expect_warning(
    loo_compare(w_list),
    "Difference in performance potentially due to chance. See McLatchie and Vehtari (2023) for details.",
    fixed = TRUE
  )

  w_list_short <- lapply(1:4, function(x) {
    suppressWarnings(waic(LLarr + rnorm(1, 0, 0.1)))
  })
  expect_no_warning(loo_compare(w_list_short))
})


comp_colnames <- c(
  "model",
  "elpd_diff",
  "se_diff",
  "p_worse",
  "diag_diff",
  "diag_elpd",
  "elpd_waic",
  "se_elpd_waic",
  "p_waic",
  "se_p_waic",
  "waic",
  "se_waic"
)

test_that("loo_compare returns expected results (2 models)", {
  comp1 <- loo_compare(w1, w1)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "data.frame")
  expect_equal(colnames(comp1), comp_colnames)
  expect_equal(comp1$model, c("model1", "model2"))
  expect_equal(comp1$elpd_diff, c(0, 0), ignore_attr = TRUE)
  expect_equal(comp1$se_diff, c(0, 0), ignore_attr = TRUE)
  expect_equal(comp1$p_worse, c(NA_real_, NA_real_), ignore_attr = TRUE)
  expect_snapshot_value(comp1, style = "serialize")
  expect_snapshot(print(comp1))

  comp2 <- loo_compare(w1, w2)
  expect_s3_class(comp2, "compare.loo")
  expect_equal(colnames(comp2), comp_colnames)
  expect_equal(comp2$p_worse, c(NA, 1))
  expect_equal(comp2$diag_diff, c("", "N < 100"))
  expect_equal(comp2$diag_elpd, c("", ""))
  expect_snapshot_value(comp2, style = "serialize")
  expect_snapshot(print(comp2))
  expect_snapshot(print(comp2, p_worse = FALSE))

  # specifying objects via ... and via arg x gives equal results
  expect_equal(comp2, loo_compare(x = list(w1, w2)))

  # custom naming works
  comp3 <- loo_compare(x = list("A" = w2, "B" = w1))
  expect_equal(comp3$model, c("B", "A"))
})


test_that("loo_compare returns expected result (3 models)", {
  w3 <- suppressWarnings(waic(LLarr3))
  comp1 <- loo_compare(w1, w2, w3)

  expect_equal(colnames(comp1), comp_colnames)
  expect_equal(comp1$model, c("model1", "model2", "model3"))
  expect_equal(comp1$p_worse, c(NA, 1, 1))
  expect_equal(comp1$diag_diff, c("", "N < 100", "N < 100"))
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "data.frame")
  expect_snapshot_value(comp1, style = "serialize")
  expect_snapshot(print(comp1))

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'x' argument
  expect_equal(comp1, loo_compare(x = list(w1, w2, w3)), ignore_attr = TRUE)
})

# Tests for deprecated compare() ------------------------------------------

test_that("compare throws deprecation warnings", {
  expect_warning(loo::compare(w1, w2), "Deprecated")
  expect_warning(loo::compare(w1, w1, w2), "Deprecated")
})

test_that("compare returns expected result (2 models)", {
  expect_warning(comp1 <- loo::compare(w1, w1), "Deprecated")
  expect_equal(comp1[1:2], c(elpd_diff = 0, se = 0))

  expect_warning(comp2 <- loo::compare(w1, w2), "Deprecated")
  expect_equal(round(comp2[1:2], 3), c(elpd_diff = -4.057, se = 0.088))
  expect_s3_class(comp2, "old_compare.loo")

  # specifying objects via ... and via arg x gives equal results
  expect_warning(comp_via_list <- loo::compare(x = list(w1, w2)), "Deprecated")
  expect_equal(comp2, comp_via_list)
})

test_that("compare returns expected result (3 models)", {
  w3 <- suppressWarnings(waic(LLarr3))
  expect_warning(comp1 <- loo::compare(w1, w2, w3), "Deprecated")

  expect_equal(
    colnames(comp1),
    c(
      "elpd_diff",
      "se_diff",
      "elpd_waic",
      "se_elpd_waic",
      "p_waic",
      "se_p_waic",
      "waic",
      "se_waic"
    )
  )
  expect_equal(rownames(comp1), c("w1", "w2", "w3"))
  expect_equal(comp1[1, 1], 0)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "matrix")
  expect_snapshot_value(comp1, style = "serialize")

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'x' argument
  expect_warning(
    comp_via_list <- loo::compare(x = list(w1, w2, w3)),
    "Deprecated"
  )
  expect_equal(comp1, comp_via_list, ignore_attr = TRUE)
})

test_that("compare throws appropriate errors", {
  expect_error(
    suppressWarnings(loo::compare(w1, w2, x = list(w1, w2))),
    "should not be specified"
  )
  expect_error(suppressWarnings(loo::compare(x = 2)), "must be a list")
  expect_error(
    suppressWarnings(loo::compare(x = list(2))),
    "should have class 'loo'"
  )
  expect_error(
    suppressWarnings(loo::compare(x = list(w1))),
    "requires at least two models"
  )

  w3 <- suppressWarnings(waic(LLarr2[,, -1]))
  expect_error(
    suppressWarnings(loo::compare(x = list(w1, w3))),
    "same number of data points"
  )
  expect_error(
    suppressWarnings(loo::compare(x = list(w1, w2, w3))),
    "same number of data points"
  )
})
