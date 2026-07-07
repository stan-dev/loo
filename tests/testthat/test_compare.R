set.seed(123)

LLarr <- example_loglik_array()
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- suppressWarnings(waic(LLarr))
w2 <- suppressWarnings(waic(LLarr2))

test_that("loo_compare throws appropriate errors", {
  w3 <- suppressWarnings(waic(LLarr[,, -1]))
  w4 <- suppressWarnings(waic(LLarr[,, -(1:2)]))

  expect_error(loo_compare(2, 3), "must be a list if not a 'loo' object")
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
  expect_true(all(is.na(comp$r2_se_diff)))
  expect_true(all(is.na(comp$mse_se_diff)))
  expect_false("r2_loo_diff" %in% colnames(comp))
  expect_false("mse_p_worse" %in% colnames(comp))
  expect_warning(
    print(comp),
    "se_diff unavailable for: r2, mse"
  )

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
  expect_equal(attr(pm1, "measure_revert_sign")$mse, FALSE)
  expect_equal(attr(pm1, "measure_revert_sign")$elpd, FALSE)
  expect_true(loo:::.measure_lower_is_better("mse_loo", loos))
  expect_false(loo:::.measure_lower_is_better("r2_loo", loos))
  expect_true(loo:::.measure_lower_is_better("mse_loo"))
  expect_false(loo:::.measure_lower_is_better("r2_loo"))
  expect_equal(
    loo:::.compare_sign_converted_measures(c("elpd_loo", "mse_loo", "r2_loo"), loos),
    c("mse")
  )

  diffs <- loo:::.measure_diffs(pm1, pm2, "elpd_loo")
  expect_length(diffs, nrow(pm1$pointwise))
  expect_equal(loo:::.measure_diffs(pm1, pm1, "elpd_loo"), rep(0, length(diffs)))
  pair_stats <- loo:::.pair_measure_stats(pm2, pm1, "elpd_loo", "sum")
  expect_equal(unname(pair_stats["se"]), loo:::se_elpd_diff(diffs))

  expect_equal(
    loo:::.measure_diff_from_estimates(pm2, pm1, "mse_loo"),
    pm1$estimates["mse_loo", "Estimate"] - pm2$estimates["mse_loo", "Estimate"]
  )
  expect_equal(
    loo:::.measure_diff_from_estimates(pm2, pm1, "r2_loo"),
    pm2$estimates["r2_loo", "Estimate"] - pm1$estimates["r2_loo", "Estimate"]
  )
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "mse_loo"), "estimates_only")
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "r2_loo"), "estimates_only")
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "elpd_loo"), "sum")
  pair_mse <- loo:::.pair_measure_stats(pm2, pm1, "mse_loo", "estimates_only")
  expect_true(is.na(pair_mse["se"]))
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
