set.seed(123)

LLarr <- example_loglik_array()
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- suppressWarnings(waic(LLarr))
w2 <- suppressWarnings(waic(LLarr2))

test_that("model_compare throws appropriate errors", {
  w3 <- suppressWarnings(waic(LLarr[,, -1]))
  w4 <- suppressWarnings(waic(LLarr[,, -(1:2)]))

  expect_error(model_compare(2, 3), "must be a list if not a 'loo' or 'pred_measure' object")
  expect_error(
    model_compare(w1, w2, x = list(w1, w2)),
    "If 'x' is a list then '...' should not be specified"
  )
  expect_error(model_compare(w1, list(1, 2, 3)), "class 'loo'")
  expect_error(model_compare(w1), "At least two models are required for comparison")
  expect_error(model_compare(x = list(w1)), "At least two models are required for comparison")
  expect_error(
    model_compare(w1, w3),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'model1' (32), 'model2' (31)",
    fixed = TRUE
  )
  expect_error(
    model_compare(w1, w2, w3),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'model1' (32), 'model2' (32), 'model3' (31)",
    fixed = TRUE
  )
  expect_error(
    model_compare(x = list("Model A" = w1, "Model B" = w2, "Model C" = w3)),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'Model A' (32), 'Model B' (32), 'Model C' (31)",
    fixed = TRUE
  )
})

test_that("model_compare dispatches loo_pred_measure inputs", {
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

  comp <- suppressMessages(model_compare(pm1, pm2))
  expect_s3_class(comp, "compare.loo")
  expect_equal(
    attr(comp, "rank_by"),
    list(kind = "default", measure = "elpd", model = NULL)
  )
  expect_true(all(c("elpd_diff", "se_diff", "p_worse", "diag_diff") %in% colnames(comp)))
  expect_true(all(c("r2_diff", "r2_se_diff", "mse_diff", "mse_se_diff") %in% colnames(comp)))
  expect_false(anyNA(comp$r2_se_diff))
  expect_false(anyNA(comp$mse_se_diff))
  expect_false("r2_loo_diff" %in% colnames(comp))
  expect_false("mse_p_worse" %in% colnames(comp))

  expect_error(
    model_compare(w1, pm1),
    "Cannot mix 'pred_measure' objects with plain 'loo' objects",
    fixed = TRUE
  )
  expect_error(
    model_compare(pm1),
    "At least two models are required for comparison",
    fixed = TRUE
  )
  expect_equal(
    attr(model_compare(w1, w2), "rank_by"),
    list(kind = "default", measure = "elpd", model = NULL)
  )
})

test_that("model_compare warns when predictive measures differ across models", {
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
    comp <- suppressMessages(model_compare(list(m1 = pm1, m2 = pm2))),
    "Omitted measures: mae \\(m2\\), mse \\(m1\\)"
  )
  expect_equal(attr(comp, "compare_measures"), c("elpd", "r2"))
  expect_false("mse_diff" %in% colnames(comp))
  expect_false("mae_diff" %in% colnames(comp))
})

test_that("model_compare works with three loo_pred_measure models", {
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

  comp <- model_compare(
    list("A" = pm1, "B" = pm2, "C" = pm3),
    rank_by = "mae"
  )
  expect_snapshot(print(comp))
  expect_equal(nrow(comp), 3L)
  expect_equal(comp$model, c("C", "B", "A"))
  expect_equal(
    attr(comp, "rank_by"),
    list(kind = "measure", measure = "mae", model = NULL)
  )
  expect_equal(attr(comp, "compare_measures"), c("elpd", "r2", "mae"))
  expect_equal(comp$mae_diff[1L], 0)
  expect_true(all(comp$mae_diff[-1L] < 0))
  # `rank_by` pins the mae-best model as the reference for *every* measure, so
  # only the reference row is zero; here C is both the mae-best and the
  # elpd-best model, so `elpd_diff` comes out negative for B and A
  expect_equal(comp$elpd_diff[1L], 0)
  expect_lt(comp$elpd_diff[comp$model == "B"], 0)
  expect_lt(comp$elpd_diff[comp$model == "A"], 0)
  expect_equal(attr(comp, "sign_converted_measures"), c("mae"))
})

test_that("model_compare informs when measure signs are converted", {
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

  expect_snapshot(comp <- model_compare(pm1, pm2))
  expect_equal(attr(comp, "sign_converted_measures"), "mse")

  pm_elpd <- loo_pred_measure(
    loo = res$loo_p_m1,
    y = res$y,
    mupred = res$mupred_m1,
    ylp = res$ylp_m1
  )
  expect_no_message(model_compare(pm_elpd, pm_elpd))
})

test_that("model_compare rank_by changes order for loo_pred_measure", {
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

  comp_elpd <- model_compare(pm1, pm2, rank_by = "elpd")
  comp_mse <- model_compare(pm1, pm2, rank_by = "mae")
  expect_equal(
    attr(comp_elpd, "rank_by"),
    list(kind = "measure", measure = "elpd", model = NULL)
  )
  expect_equal(
    attr(comp_mse, "rank_by"),
    list(kind = "measure", measure = "mae", model = NULL)
  )
  expect_equal(comp_elpd$elpd_diff[1L], 0)
  expect_equal(comp_mse$mae_diff[1L], 0)
})

test_that("without `rank_by` each measure uses its own best model as reference", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  mk <- function(m) {
    loo_pred_measure(
      loo = res[[paste0("loo_p_m", m)]],
      y = res$y,
      mupred = res[[paste0("mupred_m", m)]],
      ylp = res[[paste0("ylp_m", m)]],
      measure = c("r2", "mse", "mae")
    )
  }
  pms <- list(m1 = mk(1), m2 = mk(2), m3 = mk(3))

  comp <- suppressMessages(model_compare(pms))
  refs <- attr(comp, "compare_reference")
  expect_named(refs, c("elpd", "r2", "mse", "mae"), ignore.order = TRUE)

  # rows are still ordered by elpd, so the elpd reference is the first row
  expect_equal(refs[["elpd"]], comp$model[[1L]])
  expect_equal(comp$elpd_diff[[1L]], 0)

  for (measure in c("r2", "mse", "mae")) {
    diff_col <- comp[[paste0(measure, "_diff")]]
    # exactly one zero difference, at that measure's own best model
    expect_equal(sum(diff_col == 0), 1L)
    expect_equal(comp$model[[which(diff_col == 0)]], refs[[measure]])
    expect_true(all(diff_col <= 0))
  }

  # mse and elpd disagree here, which is the point of the per-measure reference
  expect_false(identical(refs[["mse"]], refs[["elpd"]]))

  # `rank_by` instead pins a single reference for every measure
  ranked <- suppressMessages(model_compare(pms, rank_by = "mse"))
  ranked_refs <- attr(ranked, "compare_reference")
  expect_true(all(ranked_refs == ranked$model[[1L]]))
  expect_equal(ranked$mse_diff[[1L]], 0)
})

test_that("each printed measure table is sorted best model first", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  mk <- function(m) {
    loo_pred_measure(
      loo = res[[paste0("loo_p_m", m)]],
      y = res$y,
      mupred = res[[paste0("mupred_m", m)]],
      ylp = res[[paste0("ylp_m", m)]],
      measure = c("r2", "mse", "mae")
    )
  }
  pms <- list(m1 = mk(1), m2 = mk(2), m3 = mk(3))

  printed_order <- function(comp, measure) {
    out <- utils::capture.output(
      suppressMessages(print(comp, measures = measure))
    )
    # Drop everything above the measure's own table: the PSIS-LOO diagnostics
    # block lists model names too, but says nothing about measure ordering.
    header <- grep(paste0("^-- ", measure, " "), out)
    out <- out[seq.int(header[[1L]] + 1L, length(out))]
    rows <- out[grepl("^\\s+m[0-9]", out)]
    sub("^\\s*(\\S+).*$", "\\1", rows)
  }

  for (comp in list(
    suppressMessages(model_compare(pms)),
    suppressMessages(model_compare(pms, rank_by = "mse"))
  )) {
    for (measure in c("elpd", "r2", "mse", "mae")) {
      diff_col <- if (measure == "elpd") "elpd_diff" else paste0(measure, "_diff")
      ord <- order(comp[[diff_col]], decreasing = TRUE)
      expect_equal(printed_order(comp, measure), comp$model[ord])
      # the best model on the measure leads, and the table runs downhill
      expect_equal(ord[[1L]], which.max(comp[[diff_col]]))
      expect_false(is.unsorted(rev(comp[[diff_col]][ord])))
    }
  }
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

  comp <- suppressMessages(model_compare(list(m1 = pm1, m2 = pm2, m3 = pm3)))
  expect_snapshot(print(comp))
  expect_snapshot(print(comp, measures = "all", digits = 2))
  expect_snapshot(print(comp, measures = c("r2", "mae")))

  comp_mae <- suppressMessages(model_compare(list(m1 = pm1, m2 = pm2), rank_by = "mae"))
  expect_snapshot(print(comp_mae))

  expect_error(
    print(comp, measures = "foo"),
    "Unknown measure\\(s\\) in `measures`"
  )
})

test_that("model_compare measure helpers work as expected", {
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
  expect_equal(attr(pm1, "measure_info")$elpd$diff_method, "sum")
  expect_false(attr(pm1, "measure_info")$elpd$loss)
  expect_false(attr(pm1, "measure_info")$r2$loss)
  expect_equal(attr(pm1, "measure_info")$mse$loss, TRUE)
  expect_equal(attr(pm1, "measure_info")$mse$diff_method, "mean")
  expect_equal(attr(pm1, "measure_info")$r2$diff_method, "measure_specific")
  expect_equal(attr(pm1, "measure_info")$r2$se_diff_fun, "r2")
  expect_equal(
    attr(pm1, "measure_info")$r2$extra$mse_y_i,
    (res$y - mean(res$y))^2
  )
  # only measures that need it carry `extra`
  expect_null(attr(pm1, "measure_info")$mse$extra)
  expect_true(loo:::.measure_is_loss("mse_loo", loos))
  expect_false(loo:::.measure_is_loss("r2_loo", loos))
  expect_true(loo:::.measure_is_loss("mse_loo"))
  expect_false(loo:::.measure_is_loss("r2_loo"))
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
    pm2, pm1, "r2_loo", "measure_specific", loos = loos
  )
  expect_equal(
    unname(pair_r2["diff"]),
    pm2$estimates["r2_loo", "Estimate"] - pm1$estimates["r2_loo", "Estimate"]
  )
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "mse_loo"), "mean")
  expect_equal(loo:::.measure_pointwise_diff_method(loos, "r2_loo"), "measure_specific")
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

  expect_equal(loo:::.measure_pointwise_diff_method(loos, "rmse_loo"), "measure_specific")

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

  comp <- suppressMessages(model_compare(pm1, pm2))
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

  expect_equal(loo:::.measure_pointwise_diff_method(loos, "r2_loo"), "measure_specific")

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

  comp <- suppressMessages(model_compare(pm1, pm2))
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
    meta <- attr(x, "measure_info")
    meta$r2$extra <- NULL
    attr(x, "measure_info") <- meta
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

  comp <- suppressMessages(model_compare(stale1, stale2))
  expect_false(anyNA(comp$r2_diff))
  expect_true(all(is.na(comp$r2_se_diff)))

  # one stale model does not cost the others their standard error: the
  # baseline is shared, so the other model's copy is used, and the metadata
  # check ignores `extra` rather than reporting it as a disagreement
  mixed <- suppressMessages(model_compare(stale1, pm2))
  expect_false(anyNA(mixed$r2_se_diff))
})

# two balanced-accuracy measures over the same three-class outcome: the second
# model has probability mass shifted towards the first (and largest) class, so
# the two disagree on a subset of observations and the class strata are
# unbalanced, which is where balanced accuracy differs from plain accuracy
.make_bacc_pms <- function(bias = 0.6) {
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
      measure = "bacc"
    ))
  }
  list(pm1 = make(res$mupred), pm2 = make(biased), y = y)
}

test_that("bacc differences use the stratified paired standard error", {
  fx <- .make_bacc_pms()
  pm1 <- fx$pm1
  pm2 <- fx$pm2
  loos <- list(pm1, pm2)

  expect_equal(loo:::.measure_pointwise_diff_method(loos, "bacc_loo"), "measure_specific")
  expect_equal(
    attr(pm1, "measure_info")$bacc$se_diff_fun,
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
  class_id <- attr(pm1, "measure_info")$bacc$extra$class_id
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

  comp <- suppressMessages(model_compare(pm1, pm2))
  expect_false(anyNA(comp$bacc_se_diff))
  expect_equal(comp$bacc_se_diff[1], 0)
})

test_that("bacc reports the difference without an se when the strata are gone", {
  fx <- .make_bacc_pms()

  # objects computed before the class strata were stored
  drop_strata <- function(x) {
    meta <- attr(x, "measure_info")
    meta$bacc$extra <- NULL
    attr(x, "measure_info") <- meta
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
  mixed <- suppressMessages(model_compare(fx$pm1, stale2))
  expect_false(anyNA(mixed$bacc_se_diff))

  # but a stale model at the head of the table has no second copy to fall back
  # on for its own row, which is a comparison against itself
  mixed_stale_first <- suppressMessages(model_compare(stale1, fx$pm2))
  expect_true(is.na(mixed_stale_first$bacc_se_diff[1L]))
  expect_false(is.na(mixed_stale_first$bacc_se_diff[2L]))
})

test_that("custom measures take their se_diff from `custom_se_fn`", {
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

  pms <- list(
    m1 = make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_rmse),
    m2 = make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_rmse)
  )

  # custom measures declare nothing about their standard error
  expect_equal(
    attr(pms[[1L]], "measure_info")$my_rmse$diff_method,
    "custom"
  )
  expect_null(attr(pms[[1L]], "measure_info")$my_rmse$se_diff_fun)
  expect_equal(
    loo:::.measure_pointwise_diff_method(pms, "my_rmse_loo"),
    "custom"
  )

  # omitting `custom_se_fn` is an error that names the measure
  expect_error(
    suppressMessages(model_compare(pms)),
    "my_rmse.*custom measure.*must be supplied"
  )

  # an explicit NULL reports the difference with an NA standard error
  comp_null <- suppressMessages(model_compare(pms, custom_se_fn = NULL))
  expect_false(is.na(comp_null$my_rmse_diff[[2L]]))
  expect_true(all(is.na(comp_null$my_rmse_se_diff)))
  expect_true(is.na(
    loo:::.pair_measure_stats(pms[[2L]], pms[[1L]], "my_rmse_loo", loos = pms)["se"]
  ))

  # a function gets the delta-method standard error
  comp_fn <- suppressMessages(
    model_compare(pms, custom_se_fn = loo:::.se_diff_rmse)
  )
  expect_false(any(is.na(comp_fn$my_rmse_se_diff)))
  # without `rank_by`, `my_rmse` is compared against its own best model, which
  # is the row with a zero difference and a zero standard error
  ref_name <- attr(comp_fn, "compare_reference")[["my_rmse"]]
  cmp_name <- setdiff(names(pms), ref_name)
  ref_row <- match(ref_name, comp_fn$model)
  cmp_row <- match(cmp_name, comp_fn$model)
  expect_equal(comp_fn$my_rmse_se_diff[[ref_row]], 0)
  expect_gt(comp_fn$my_rmse_se_diff[[cmp_row]], 0)
  # the difference itself does not depend on how the SE was obtained
  expect_equal(comp_fn$my_rmse_diff, comp_null$my_rmse_diff)

  pair_ref <- loo:::.pair_measure_stats(
    pms[[cmp_name]], pms[[ref_name]], "my_rmse_loo",
    loos = pms, se_fn = loo:::.se_diff_rmse
  )
  expect_equal(unname(pair_ref["se"]), comp_fn$my_rmse_se_diff[[cmp_row]])

  pair <- loo:::.pair_measure_stats(
    pms[[2L]], pms[[1L]], "my_rmse_loo",
    loos = pms, se_fn = loo:::.se_diff_rmse
  )

  # a custom measure can carry its own auxiliary data through to `custom_se_fn`,
  # and each model receives its own copy
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
  scaled_se_fn <- function(ref, cmp) {
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
    attr(scaled[[1L]], "measure_info")$my_scaled$extra,
    list(scale = stats::sd(res$y), n_used = length(res$y))
  )
  pair_scaled <- loo:::.pair_measure_stats(
    scaled[[2L]], scaled[[1L]], "my_scaled_loo",
    loos = scaled, se_fn = scaled_se_fn
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

  # a `custom_se_fn` that returns nonsense is caught
  expect_error(
    suppressMessages(
      model_compare(pms, custom_se_fn = function(ref, cmp) c(1, 2))
    ),
    "must return a numeric scalar"
  )
})

test_that("a declared custom loss is compared and ranked as a loss", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")

  # squared error, whose estimate is the mean of its pointwise values
  make_fun <- function(declare_loss) {
    f <- function(y, mupred) {
      sqe <- (y - colMeans(mupred))^2
      list(
        estimate = mean(sqe),
        se = sqrt(var(sqe) / length(sqe)),
        pointwise = sqe
      )
    }
    attr(f, "measure_name") <- "my_mse"
    if (declare_loss) attr(f, "measure_loss") <- TRUE
    f
  }

  make <- function(m, fun) {
    loo_pred_measure(
      loo = res[[paste0("loo_p_m", m)]],
      y = res$y,
      mupred = res[[paste0("mupred_m", m)]],
      ylp = res[[paste0("ylp_m", m)]],
      measure = fun
    )
  }

  declared <- list(m1 = make(1, make_fun(TRUE)), m2 = make(2, make_fun(TRUE)))
  plain <- list(m1 = make(1, make_fun(FALSE)), m2 = make(2, make_fun(FALSE)))

  expect_true(attr(declared$m1, "measure_info")$my_mse$loss)
  expect_true(loo:::.measure_is_loss("my_mse_loo", declared))
  expect_false(loo:::.measure_is_loss("my_mse_loo", plain))

  # the sign conversion is announced, as it is for built-in loss measures
  expect_message(
    comp <- model_compare(declared, custom_se_fn = "mean"),
    "my_mse.*utility scale"
  )
  comp_plain <- suppressMessages(model_compare(plain, custom_se_fn = "mean"))

  expect_equal(attr(comp, "sign_converted_measures"), "my_mse")
  expect_length(attr(comp_plain, "sign_converted_measures"), 0L)

  # same models, same measure: the declared loss and the undeclared utility
  # disagree about which model is best, so each picks the other's reference
  expect_false(identical(
    attr(comp, "compare_reference")[["my_mse"]],
    attr(comp_plain, "compare_reference")[["my_mse"]]
  ))
  # against a single pinned reference, only the orientation of the difference
  # changes
  comp_ref <- suppressMessages(
    model_compare(declared, rank_by = "elpd", custom_se_fn = "mean")
  )
  comp_plain_ref <- suppressMessages(
    model_compare(plain, rank_by = "elpd", custom_se_fn = "mean")
  )
  expect_equal(comp_ref$my_mse_diff, -comp_plain_ref$my_mse_diff)
  expect_equal(comp_ref$my_mse_se_diff, comp_plain_ref$my_mse_se_diff)
  # ... and the declared version agrees with the built-in `mse` on which model
  # is worse
  builtin <- list(m1 = make(1, "mse"), m2 = make(2, "mse"))
  comp_builtin <- suppressMessages(model_compare(builtin))
  expect_equal(comp$model, comp_builtin$model)
  expect_equal(sign(comp$my_mse_diff), sign(comp_builtin$mse_diff))

  # `rank_by` puts the lowest loss first
  ranked <- suppressMessages(
    model_compare(declared, rank_by = "my_mse", custom_se_fn = NULL)
  )
  ranked_plain <- suppressMessages(
    model_compare(plain, rank_by = "my_mse", custom_se_fn = NULL)
  )
  est <- vapply(
    ranked$model,
    function(m) declared[[m]]$estimates["my_mse_loo", "Estimate"],
    numeric(1)
  )
  expect_false(is.unsorted(est))
  expect_equal(rev(ranked$model), ranked_plain$model)

  # models must agree on the declaration
  expect_error(
    suppressMessages(
      model_compare(list(declared$m1, plain$m2), custom_se_fn = "mean")
    ),
    "disagree on `measure_info`"
  )
})

test_that("`custom_se_fn` accepts the \"sum\" and \"mean\" shorthands", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")

  # a custom measure reproducing the built-in `mae` on the utility scale. It
  # declares `log_weights` so that it gets the same PSIS-weighted point
  # predictions the built-in uses, and negates so that it is a genuine utility
  # (a custom measure that does not declare `measure_loss` is one).
  my_mae <- function(y, mupred, log_weights) {
    w <- exp(loo:::.normalize_and_validate_log_weights(
      log_weights = log_weights,
      n_draws = nrow(mupred),
      n_obs = ncol(mupred)
    ))
    ae_i <- -abs(y - colSums(w * mupred))
    list(
      estimate = mean(ae_i),
      se = sqrt(var(ae_i) / length(ae_i)),
      pointwise = ae_i
    )
  }
  attr(my_mae, "measure_name") <- "my_mae"

  make <- function(loo, mupred, ylp, measure) {
    loo_pred_measure(
      loo = loo, y = res$y, mupred = mupred, ylp = ylp, measure = measure
    )
  }

  custom <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_mae),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_mae)
  )
  builtin <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, "mae"),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, "mae")
  )

  comp_custom <- suppressMessages(
    model_compare(custom, custom_se_fn = "mean")
  )
  comp_builtin <- suppressMessages(model_compare(builtin))

  # "mean" reuses the built-in branch, so results must match `mae` exactly
  expect_equal(comp_custom$my_mae_diff, comp_builtin$mae_diff)
  expect_equal(comp_custom$my_mae_se_diff, comp_builtin$mae_se_diff)

  # "sum" against a custom measure whose estimate is a sum of pointwise values
  my_sum <- function(y, mupred) {
    ae_i <- -abs(y - colMeans(mupred))
    list(estimate = sum(ae_i), se = sqrt(length(ae_i) * var(ae_i)),
         pointwise = ae_i)
  }
  attr(my_sum, "measure_name") <- "my_sum"
  # named so that the comparison's row order can be mapped back to the inputs
  summed <- list(
    a = make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_sum),
    b = make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_sum)
  )
  comp_sum <- suppressMessages(model_compare(summed, custom_se_fn = "sum"))
  expect_false(any(is.na(comp_sum$my_sum_se_diff)))

  # matches `sqrt(N) * sd(d_i)` computed by hand from the pointwise columns
  ref_pw <- summed[[comp_sum$model[[1L]]]]$pointwise[, "my_sum_loo"]
  cmp_pw <- summed[[comp_sum$model[[2L]]]]$pointwise[, "my_sum_loo"]
  d <- cmp_pw - ref_pw
  expect_equal(comp_sum$my_sum_se_diff[[2L]], sqrt(length(d)) * sd(d))
  expect_equal(comp_sum$my_sum_diff[[2L]], sum(d))

  # a declared aggregation that does not reproduce the estimate warns
  my_rmse <- function(y, mupred) {
    sqe_i <- (y - colMeans(mupred))^2
    list(
      estimate = sqrt(mean(sqe_i)),
      se = sqrt(var(sqe_i) / length(sqe_i)) / (2 * sqrt(mean(sqe_i))),
      pointwise = sqe_i
    )
  }
  attr(my_rmse, "measure_name") <- "my_rmse"
  mismatched <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_rmse),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_rmse)
  )
  expect_warning(
    suppressMessages(model_compare(mismatched, custom_se_fn = "mean")),
    "does not reproduce its estimate"
  )

  # any other string is rejected
  expect_error(
    suppressMessages(model_compare(custom, custom_se_fn = "median")),
    "must be a function"
  )
})

test_that("`custom_se_fn` validates its per-measure form", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")

  make_fun <- function(name, offset) {
    f <- function(y, mupred) {
      ae_i <- abs(y - colMeans(mupred)) + offset
      list(estimate = mean(ae_i), se = sqrt(var(ae_i) / length(ae_i)),
           pointwise = ae_i)
    }
    attr(f, "measure_name") <- name
    f
  }
  a <- make_fun("m_a", 0)
  b <- make_fun("m_b", 1)

  make <- function(loo, mupred, ylp, measure) {
    loo_pred_measure(
      loo = loo, y = res$y, mupred = mupred, ylp = ylp, measure = measure
    )
  }

  two <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, list(m_a = a, m_b = b)),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, list(m_a = a, m_b = b))
  )

  # a bare value is ambiguous with more than one custom measure
  expect_error(
    suppressMessages(model_compare(two, custom_se_fn = "mean")),
    "must be a named list"
  )

  # a named list may mix the accepted forms
  comp <- suppressMessages(model_compare(
    two,
    custom_se_fn = list(m_a = "mean", m_b = NULL)
  ))
  expect_false(any(is.na(comp$m_a_se_diff)))
  expect_true(all(is.na(comp$m_b_se_diff)))

  # an entry must exist for every custom measure
  expect_error(
    suppressMessages(model_compare(two, custom_se_fn = list(m_a = "mean"))),
    "no entry for custom measure"
  )
  # unknown names are typos
  expect_error(
    suppressMessages(model_compare(
      two,
      custom_se_fn = list(m_a = "mean", m_b = NULL, nope = "mean")
    )),
    "Unknown measure"
  )
  # unnamed lists cannot be matched to measures
  expect_error(
    suppressMessages(model_compare(two, custom_se_fn = list("mean", NULL))),
    "must be named after a custom measure"
  )
  # elements must be one of the accepted forms
  expect_error(
    suppressMessages(model_compare(
      two,
      custom_se_fn = list(m_a = 1, m_b = NULL)
    )),
    "must be a function"
  )

  # supplying it when no custom measure is present warns and changes nothing
  builtin <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, "mae"),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, "mae")
  )
  expect_warning(
    comp_builtin <- suppressMessages(
      model_compare(builtin, custom_se_fn = "mean")
    ),
    "only used for custom measures"
  )
  expect_equal(
    comp_builtin$mae_se_diff,
    suppressMessages(model_compare(builtin))$mae_se_diff
  )

  # the model_compare() alias forwards the argument
  one <- list(
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, a),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, a)
  )
  expect_equal(
    suppressMessages(model_compare(one, custom_se_fn = "mean"))$m_a_se_diff,
    suppressMessages(model_compare(one, custom_se_fn = "mean"))$m_a_se_diff
  )
  expect_error(
    suppressMessages(model_compare(one)),
    "must be supplied"
  )
})

test_that("model_compare errors on inconsistent measure metadata", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  # the same custom measure, but only one model declares it a loss
  make_fun <- function(loss) {
    fun <- function(y, mupred) {
      e <- (y - colMeans(mupred))^2
      list(estimate = mean(e), se = sd(e) / sqrt(length(e)), pointwise = e)
    }
    attr(fun, "measure_loss") <- loss
    fun
  }
  make <- function(loo, mupred, ylp, fun) {
    loo_pred_measure(
      loo = loo,
      y = res$y,
      mupred = mupred,
      ylp = ylp,
      measure = list(my_mse = fun)
    )
  }
  pm1 <- make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, make_fun(TRUE))
  pm2 <- make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, make_fun(FALSE))

  expect_error(
    suppressMessages(model_compare(pm1, pm2, custom_se_fn = list(my_mse = "mean"))),
    "disagree on `measure_info` for measure 'my_mse'"
  )
})

test_that("model_compare errors when compare metadata is missing on some models", {
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
  measure_info <- attr(pm2, "measure_info")
  measure_info$mse <- NULL
  attr(pm2, "measure_info") <- measure_info

  expect_error(
    suppressMessages(model_compare(pm1, pm2)),
    "Not all models provide `measure_info` for measure 'mse'"
  )
})

test_that("model_compare warns when rank_by is ignored for classic loo objects", {
  expect_warning(
    model_compare(w1, w2, rank_by = "mse"),
    "`rank_by` is only used for `pred_measure` comparisons"
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

test_that("model_compare works for all built-in measures", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  res_roaches <- readRDS("data-for-tests/test_data_roaches.Rds")
  roaches_measures <- c(
    "ic", "mlpd", "mae", "r2", "rmse", "mse"
  )
  for (measure in roaches_measures) {
    pm1 <- .make_compare_pm(res, 1L, measure)
    pm2 <- .make_compare_pm(res, 2L, measure)
    comp <- suppressMessages(model_compare(pm1, pm2))
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
    comp <- suppressMessages(model_compare(pm1, pm2))
    expect_true(paste0(measure, "_diff") %in% colnames(comp), info = measure)
    expect_equal(attr(comp, "compare_measures"), c("elpd", measure), info = measure)
  }

  for (measure in c("brier", "acc", "bacc")) {
    pm1 <- .make_compare_pm_synthetic(measure)
    pm2 <- .make_compare_pm_synthetic(measure)
    comp <- suppressMessages(model_compare(pm1, pm2))
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

test_that("model_compare warns for many loo_pred_measure models", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  set.seed(123)
  pm_list <- .make_many_compare_pms(res, 25L)
  expect_warning(
    suppressMessages(model_compare(pm_list)),
    "Difference in performance potentially due to chance. See McLatchie and Vehtari (2023) for details.",
    fixed = TRUE
  )

  pm_list_short <- .make_many_compare_pms(res, 4L)
  expect_no_warning(suppressMessages(model_compare(pm_list_short)))
})

test_that("model_compare throws appropriate warnings", {
  w3 <- w1
  w4 <- w2
  class(w3) <- class(w4) <- c("kfold", "loo")
  attr(w3, "K") <- 2
  attr(w4, "K") <- 3
  expect_warning(
    model_compare(w3, w4),
    "Not all kfold objects have the same K value"
  )

  class(w4) <- c("psis_loo", "loo")
  attr(w4, "K") <- NULL
  expect_warning(model_compare(w3, w4), "Comparing LOO-CV to K-fold-CV")

  w3 <- w1
  w4 <- w2
  attr(w3, "yhash") <- "a"
  attr(w4, "yhash") <- "b"
  expect_warning(model_compare(w3, w4), "Not all models have the same y variable")

  set.seed(123)
  w_list <- lapply(1:25, function(x) {
    suppressWarnings(waic(LLarr + rnorm(1, 0, 0.1)))
  })
  expect_warning(
    model_compare(w_list),
    "Difference in performance potentially due to chance. See McLatchie and Vehtari (2023) for details.",
    fixed = TRUE
  )

  w_list_short <- lapply(1:4, function(x) {
    suppressWarnings(waic(LLarr + rnorm(1, 0, 0.1)))
  })
  expect_no_warning(model_compare(w_list_short))
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

test_that("model_compare returns expected results (2 models)", {
  comp1 <- model_compare(w1, w1)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "data.frame")
  expect_equal(colnames(comp1), comp_colnames)
  expect_equal(comp1$model, c("model1", "model2"))
  expect_equal(comp1$elpd_diff, c(0, 0), ignore_attr = TRUE)
  expect_equal(comp1$se_diff, c(0, 0), ignore_attr = TRUE)
  expect_equal(comp1$p_worse, c(NA_real_, NA_real_), ignore_attr = TRUE)
  expect_snapshot_value(comp1, style = "serialize")
  expect_snapshot(print(comp1))

  comp2 <- model_compare(w1, w2)
  expect_s3_class(comp2, "compare.loo")
  expect_equal(colnames(comp2), comp_colnames)
  expect_equal(comp2$p_worse, c(NA, 1))
  expect_equal(comp2$diag_diff, c("", "N < 100"))
  expect_equal(comp2$diag_elpd, c("", ""))
  expect_snapshot_value(comp2, style = "serialize")
  expect_snapshot(print(comp2))
  expect_snapshot(print(comp2, p_worse = FALSE))
  expect_snapshot(print(comp2, simplify = FALSE))
  expect_snapshot(print(comp2, simplify = FALSE, p_worse = FALSE))

  # specifying objects via ... and via arg x gives equal results
  expect_equal(comp2, model_compare(x = list(w1, w2)))

  # custom naming works
  comp3 <- model_compare(x = list("A" = w2, "B" = w1))
  expect_equal(comp3$model, c("B", "A"))
})

test_that("model_compare returns expected result (3 models)", {
  w3 <- suppressWarnings(waic(LLarr3))
  comp1 <- model_compare(w1, w2, w3)

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
  expect_equal(comp1, model_compare(x = list(w1, w2, w3)), ignore_attr = TRUE)
})

test_that("model_compare with simplify=FALSE returns expected result", {
  LL <- example_loglik_array()
  loo1 <- loo(LL)
  loo2 <- loo(LL + 1)
  loo3 <- loo(LL + 2)
  comp <- model_compare(loo1, loo2, loo3)
  expect_snapshot(print(comp, simplify = FALSE))
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

# model_compare across evaluation sources -----------------------------------

.compare_src_res <- function() readRDS("data-for-tests/test_data_roaches.Rds")

.jitter_mupred <- function(mupred, sd) {
  mupred + stats::rnorm(length(mupred), 0, sd)
}

test_that("model_compare compares kfold_pred_measure objects", {
  res <- .compare_src_res()
  set.seed(4321)
  k1 <- kfold_pred_measure(
    y = res$y, mupred = res$mupred, kfold = res$kfold,
    measure = c("rmse", "mse")
  )
  k2 <- kfold_pred_measure(
    y = res$y, mupred = .jitter_mupred(res$mupred, 3), kfold = res$kfold,
    measure = c("rmse", "mse")
  )

  comp <- suppressMessages(model_compare(list(m1 = k1, m2 = k2)))
  expect_s3_class(comp, "compare.loo")
  expect_equal(attr(comp, "compare_source"), "kfold")
  expect_equal(attr(comp, "compare_measures"), c("elpd", "rmse", "mse"))

  # measures are matched on bare names, with the `_kfold` suffix stripped
  expect_true(all(c("rmse_diff", "rmse_se_diff", "mse_diff", "mse_se_diff") %in%
                    colnames(comp)))
  expect_false(any(grepl("_kfold_diff$", colnames(comp))))
  expect_false(anyNA(comp$rmse_se_diff))

  # Pareto k diagnostics do not exist outside PSIS-LOO
  expect_false("diag_elpd" %in% colnames(comp))
  # ELPD-only columns are still produced
  expect_true(all(c("elpd_diff", "se_diff", "p_worse", "diag_diff") %in%
                    colnames(comp)))
})

test_that("model_compare compares test_pred_measure objects", {
  res_cv <- readRDS("data-for-tests/test_data_sleep_cv.Rds")
  set.seed(4321)
  t1 <- test_pred_measure(
    y = res_cv$y_test, mupred = res_cv$mupred_test,
    ylp_test = res_cv$ylp_test, measure = "rmse"
  )
  t2 <- test_pred_measure(
    y = res_cv$y_test, mupred = .jitter_mupred(res_cv$mupred_test, 5),
    ylp_test = res_cv$ylp_test, measure = "rmse"
  )

  comp <- suppressMessages(model_compare(list(m1 = t1, m2 = t2)))
  expect_equal(attr(comp, "compare_source"), "test")
  expect_equal(attr(comp, "compare_measures"), c("elpd", "rmse"))
  expect_true(all(c("rmse_diff", "rmse_se_diff") %in% colnames(comp)))
  expect_false("diag_elpd" %in% colnames(comp))
})

test_that("model_compare warns that insample_pred_measure comparisons are biased", {
  res <- .compare_src_res()
  set.seed(4321)
  i1 <- insample_pred_measure(
    y = res$y, mupred = res$mupred, ylp = res$ylp, measure = "rmse"
  )
  i2 <- insample_pred_measure(
    y = res$y, mupred = .jitter_mupred(res$mupred, 3), ylp = res$ylp,
    measure = "rmse"
  )

  expect_warning(
    comp <- suppressMessages(model_compare(list(m1 = i1, m2 = i2))),
    "optimistically biased"
  )
  expect_equal(attr(comp, "compare_source"), "insample")
  # in-sample measures carry no suffix at all
  expect_equal(attr(comp, "compare_measures"), c("elpd", "rmse"))
  expect_true(all(c("rmse_diff", "rmse_se_diff") %in% colnames(comp)))
})

test_that("model_compare errors when evaluation sources are mixed", {
  res <- .compare_src_res()
  k1 <- kfold_pred_measure(y = res$y, mupred = res$mupred, kfold = res$kfold,
                           measure = "rmse")
  l1 <- loo_pred_measure(loo = res$loo, y = res$y, mupred = res$mupred,
                         measure = "rmse")
  i1 <- insample_pred_measure(y = res$y, mupred = res$mupred, ylp = res$ylp,
                              measure = "rmse")

  # all three have the same number of observations, so this is genuinely the
  # source check firing rather than the observation-count check
  expect_equal(nrow(k1$pointwise), nrow(l1$pointwise))
  expect_error(
    model_compare(k1, l1),
    "All models must be evaluated on the same source",
    fixed = TRUE
  )
  expect_error(
    model_compare(l1, i1),
    "All models must be evaluated on the same source",
    fixed = TRUE
  )
})

test_that("model_compare warns when kfold results use different K", {
  res <- .compare_src_res()
  set.seed(4321)
  k1 <- kfold_pred_measure(y = res$y, mupred = res$mupred, kfold = res$kfold,
                           measure = "rmse")
  k2 <- kfold_pred_measure(y = res$y, mupred = .jitter_mupred(res$mupred, 3),
                           kfold = res$kfold, measure = "rmse")
  attr(k2, "K") <- 5L

  expect_warning(
    suppressMessages(model_compare(list(m1 = k1, m2 = k2))),
    "Not all kfold objects have the same K value"
  )
})

test_that("model_compare rank_by resolves bare names for suffixed measures", {
  res <- .compare_src_res()
  set.seed(4321)
  k1 <- kfold_pred_measure(y = res$y, mupred = res$mupred, kfold = res$kfold,
                           measure = c("rmse", "mae"))
  k2 <- kfold_pred_measure(y = res$y, mupred = .jitter_mupred(res$mupred, 3),
                           kfold = res$kfold, measure = c("rmse", "mae"))

  comp <- suppressMessages(model_compare(list(m1 = k1, m2 = k2), rank_by = "mae"))
  expect_equal(
    attr(comp, "rank_by"),
    list(kind = "measure", measure = "mae", model = NULL)
  )
  expect_equal(comp$mae_diff[1L], 0)
  expect_true(all(comp$mae_diff[-1L] <= 0))

  expect_error(
    suppressMessages(model_compare(list(m1 = k1, m2 = k2), rank_by = "nope")),
    "`rank_by` value 'nope' is neither a measure nor a model name",
    fixed = TRUE
  )
})

test_that("model_compare rank_by accepts a model name as the reference model", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  mk <- function(m) {
    loo_pred_measure(
      loo = res[[paste0("loo_p_m", m)]],
      y = res$y,
      mupred = res[[paste0("mupred_m", m)]],
      ylp = res[[paste0("ylp_m", m)]],
      measure = c("r2", "mse", "mae")
    )
  }
  pms <- list(m1 = mk(1), m2 = mk(2), m3 = mk(3))

  default <- suppressMessages(model_compare(pms))
  pinned <- suppressMessages(model_compare(pms, rank_by = "m1"))

  # the named model is the reference for every measure, whether or not it is
  # the best model
  expect_equal(
    attr(pinned, "rank_by"),
    list(kind = "model", measure = "elpd", model = "m1")
  )
  expect_true(all(attr(pinned, "compare_reference") == "m1"))
  for (col in c("elpd_diff", "r2_diff", "mse_diff", "mae_diff")) {
    expect_equal(pinned[[col]][pinned$model == "m1"], 0)
  }

  # rows stay ordered by elpd, as without `rank_by`
  expect_equal(pinned$model, default$model)

  # differences are the same comparisons, just re-referenced
  expect_equal(
    pinned$elpd_diff - pinned$elpd_diff[pinned$model == default$model[[1L]]],
    default$elpd_diff
  )

  expect_output(print(pinned), "All measures compared against model m1")

  # `diag_diff` flags a *small* difference, so a large positive one --- which
  # only arises when the reference is not the best model --- stays unflagged
  large_positive <- pinned$elpd_diff[pinned$elpd_diff > 4]
  expect_true(length(large_positive) > 0)
  expect_equal(pinned$diag_diff[pinned$elpd_diff > 4], rep("", length(large_positive)))
})

test_that("diag_diff flags the magnitude of elpd_diff, not its sign", {
  expect_equal(diag_diff(500, c(0, -2, 2, -10, 10)),
               c("", "|elpd_diff| < 4", "|elpd_diff| < 4", "", ""))
  # small N takes priority over the difference itself, for every non-reference
  expect_equal(diag_diff(50, c(0, -2, 10)), c("", "N < 100", "N < 100"))
})

test_that("printed comparison output stays within 80 columns", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  mk <- function(m) {
    loo_pred_measure(
      loo = res[[paste0("loo_p_m", m)]],
      y = res$y,
      mupred = res[[paste0("mupred_m", m)]],
      ylp = res[[paste0("ylp_m", m)]],
      measure = c("r2", "mse", "mae")
    )
  }
  # long names stress both the wrapped sentences and the table layout
  nms <- c(
    "poisson_baseline_model",
    "negbin_pretreatment_model",
    "poisson_full_interaction_model"
  )
  pms <- stats::setNames(lapply(1:3, mk), nms)
  comp <- suppressMessages(model_compare(pms))

  for (measures in list(NULL, "all", c("r2", "mae"))) {
    out <- utils::capture.output(
      suppressMessages(print(comp, measures = measures))
    )
    expect_true(all(nchar(out) <= 80))
  }
})

test_that("model_compare rank_by model name works for plain loo objects", {
  comp <- model_compare(list(a = w1, b = w2), rank_by = "b")
  expect_equal(
    attr(comp, "rank_by"),
    list(kind = "model", measure = "elpd", model = "b")
  )
  expect_equal(attr(comp, "compare_reference"), c(elpd = "b"))
  expect_equal(comp$elpd_diff[comp$model == "b"], 0)
  expect_true(is.na(comp$p_worse[comp$model == "b"]))

  default <- model_compare(list(a = w1, b = w2))
  expect_equal(comp$model, default$model)
  expect_equal(
    comp$elpd_diff - comp$elpd_diff[comp$model == default$model[[1L]]],
    default$elpd_diff
  )
  expect_message(print(comp), "Differences computed against model b")
})

test_that("model_compare rank_by prefers the measure when a model shares its name", {
  res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
  mk <- function(m) {
    loo_pred_measure(
      loo = res[[paste0("loo_p_m", m)]],
      y = res$y,
      mupred = res[[paste0("mupred_m", m)]],
      ylp = res[[paste0("ylp_m", m)]],
      measure = c("mse")
    )
  }
  pms <- list(mse = mk(1), m2 = mk(2))

  expect_warning(
    comp <- suppressMessages(model_compare(pms, rank_by = "mse")),
    "matches both a measure and a model name"
  )
  expect_equal(
    attr(comp, "rank_by"),
    list(kind = "measure", measure = "mse", model = NULL)
  )

  expect_error(
    suppressMessages(model_compare(pms, rank_by = 1)),
    "`rank_by` must be a single measure name or model name",
    fixed = TRUE
  )
})

# Tests for deprecated loo_compare() --------------------------------------

test_that("loo_compare throws a deprecation warning once per session", {
  # forget that the warning was already issued earlier in this session
  forget_warning <- function() {
    assign("loo_compare", FALSE, envir = environment(loo:::.deprecate_once)$state)
  }
  forget_warning()
  on.exit(forget_warning(), add = TRUE)

  expect_warning(loo_compare(w1, w2), "deprecated")
  # already warned in this session, so these are silent
  expect_no_warning(loo_compare(w1, w2))
  expect_no_warning(loo_compare(x = list(w1, w2)))

  forget_warning()
  expect_warning(loo_compare(x = list(w1, w2)), "deprecated")
})

test_that("loo_compare still returns what model_compare returns", {
  expect_identical(suppressWarnings(loo_compare(w1, w2)), model_compare(w1, w2))
  expect_identical(
    suppressWarnings(loo_compare(x = list("A" = w1, "B" = w2))),
    model_compare(x = list("A" = w1, "B" = w2))
  )
})

test_that("loo_compare is frozen to classic elpd comparison", {
  res <- .compare_src_res()
  set.seed(4321)
  k1 <- kfold_pred_measure(y = res$y, mupred = res$mupred, kfold = res$kfold,
                           measure = "rmse")
  k2 <- kfold_pred_measure(y = res$y, mupred = .jitter_mupred(res$mupred, 3),
                           kfold = res$kfold, measure = "rmse")

  expect_error(
    suppressWarnings(loo_compare(k1, k2)),
    "Use `model_compare()` to compare 'pred_measure' results",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(loo_compare(w1, w2, rank_by = "model1")),
    "`rank_by` is not supported by the deprecated `loo_compare()`",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(loo_compare(w1, w2, custom_se_fn = "mean")),
    "`custom_se_fn` is not supported by the deprecated `loo_compare()`",
    fixed = TRUE
  )
})

test_that("loo_compare is still a generic", {
  # methods registered elsewhere (e.g. brms, rstanarm) keep dispatching
  assign("loo_compare.fake_fit", function(x, ...) "dispatched", envir = globalenv())
  on.exit(rm("loo_compare.fake_fit", envir = globalenv()), add = TRUE)
  expect_identical(loo_compare(structure(list(), class = "fake_fit")), "dispatched")
})

test_that("print.compare.loo names the source for non-loo comparisons", {
  res <- .compare_src_res()
  set.seed(4321)
  k1 <- kfold_pred_measure(y = res$y, mupred = res$mupred, kfold = res$kfold,
                           measure = "rmse")
  k2 <- kfold_pred_measure(y = res$y, mupred = .jitter_mupred(res$mupred, 3),
                           kfold = res$kfold, measure = "rmse")
  comp <- suppressMessages(model_compare(list(m1 = k1, m2 = k2)))

  expect_output(print(comp), "K-fold cross-validation", fixed = TRUE)

  # LOO is the default and stays unlabelled
  l1 <- loo_pred_measure(loo = res$loo, y = res$y, mupred = res$mupred,
                         measure = "rmse")
  l2 <- loo_pred_measure(loo = res$loo, y = res$y,
                         mupred = .jitter_mupred(res$mupred, 3),
                         measure = "rmse")
  comp_loo <- suppressMessages(model_compare(list(m1 = l1, m2 = l2)))
  expect_no_match(
    paste(capture.output(print(comp_loo)), collapse = "\n"),
    "evaluated on",
    fixed = TRUE
  )
})

test_that("rps is sign-converted for comparison but srps is not", {
  set.seed(20250826)
  S <- 400L
  n <- 60L
  y <- rnorm(n)
  # the second model is the misspecified one under both scores
  good <- matrix(rnorm(S * n), nrow = S)
  bad <- matrix(rnorm(S * n, mean = 2, sd = 3), nrow = S)
  make <- function(ypred) {
    insample_pred_measure(
      y = y,
      ypred = ypred,
      ylp = matrix(dnorm(rep(y, each = S), log = TRUE), nrow = S),
      measure = c("rps", "srps")
    )
  }
  pm1 <- make(good)
  pm2 <- make(bad)
  expect_false(anyNA(pm1$estimates))
  expect_false(anyNA(pm2$estimates))

  comp <- suppressMessages(model_compare(list(m1 = pm1, m2 = pm2)))

  # the unscaled score is a loss, so it is flipped onto the utility scale; the
  # scaled score already is a utility
  expect_equal(attr(comp, "sign_converted_measures"), "rps")

  # both differences are then on a utility scale: the best model has 0 and the
  # other a non-positive difference
  expect_true(all(comp$rps_diff <= 0))
  expect_true(all(comp$srps_diff <= 0))

  # `rank_by` follows the same orientation: the well-specified model must rank
  # first under both scores
  expect_equal(
    suppressMessages(
      model_compare(list(m1 = pm1, m2 = pm2), rank_by = "rps")
    )$model[1L],
    "m1"
  )
  expect_equal(
    suppressMessages(
      model_compare(list(m1 = pm1, m2 = pm2), rank_by = "srps")
    )$model[1L],
    "m1"
  )
})
