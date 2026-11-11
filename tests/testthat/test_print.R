# load data
temp <- readRDS("data-for-tests/test_data_roaches.Rds")

measure_specs <- function(temp) {
  list(
    list(measure = "r2", args = list(y = temp$y, mupred = temp$mupred)),
    list(measure = "rmse", args = list(y = temp$y, mupred = temp$mupred)),
    list(measure = "mse", args = list(y = temp$y, mupred = temp$mupred)),
    list(measure = "mae", args = list(y = temp$y, mupred = temp$mupred)),
    list(measure = "rps", args = list(y = temp$y, ypred = temp$ypred)),
    list(measure = "srps", args = list(y = temp$y, ypred = temp$ypred)),
    list(measure = "mlpd", args = list(ylp = temp$ylp))
  )
}

run_measure_snapshots <- function(loo_start, measures) {
  loo_iter <- loo_start
  for (i in seq_along(measures)) {
    measure <- measures[[i]]
    loo_prev <- loo_iter
    call_args <- c(
      measure$args,
      list(
        measures = measure$measure,
        save_psis = TRUE
      )
    )
    if (!is.null(measure$control)) {
      call_args$control <- measure$control
    }
    if (i == 1L) {
      call_args$loo <- loo_iter
      loo_iter <- do.call(loo_pred_measure, call_args)
    } else {
      call_args$predperf <- loo_iter
      loo_iter <- do.call(pred_measure, call_args)
    }
    common_rows <- intersect(
      rownames(loo_prev$estimates),
      rownames(loo_iter$estimates)
    )
    common_cols <- intersect(
      colnames(loo_prev$pointwise),
      colnames(loo_iter$pointwise)
    )
    expect_equal(
      loo_iter$estimates[common_rows, , drop = FALSE],
      loo_prev$estimates[common_rows, , drop = FALSE],
      info = measure$measure
    )
    expect_equal(
      loo_iter$pointwise[, common_cols, drop = FALSE],
      loo_prev$pointwise[, common_cols, drop = FALSE],
      info = measure$measure
    )
    expect_snapshot_output(print(loo_iter))
  }
  loo_iter
}

test_that("loo_pred_measure print snapshots", {
  loo_ordered <- run_measure_snapshots(temp$loo, measure_specs(temp))
  loo_shuffled <- run_measure_snapshots(
    temp$loo,
    with(set.seed(0), sample(measure_specs(temp)))
  )
  expect_setequal(
    rownames(loo_ordered$estimates),
    rownames(loo_shuffled$estimates)
  )
  expect_equal(
    loo_ordered$estimates[rownames(loo_shuffled$estimates), , drop = FALSE],
    loo_shuffled$estimates
  )
})

test_that("loo_pred_measure print output with elpd", {
  x <- loo_pred_measure(
    loo = temp$loo,
    y = temp$y,
    mupred = temp$mupred,
    measures = c("elpd", "r2")
  )
  expect_snapshot_output(print(x))
})

test_that("test_pred_measure print output", {
  res <- readRDS("data-for-tests/test_data_sleep_cv.Rds")
  x <- test_pred_measure(
    y = res$y_test,
    ypred = res$ypred_test,
    mupred = res$mupred_test,
    ylp_test = res$ylp_test,
    measures = c("rmse", "r2")
  )
  expect_s3_class(x, "test_pred_measure")
  expect_snapshot_output(print(x))
})

test_that(".se_digits takes the places from the standard error", {
  expect_equal(loo:::.se_digits(0.0003), 4)
  expect_equal(loo:::.se_digits(0.045), 3)
  expect_equal(loo:::.se_digits(1.4), 1)
  # a large SE still gets `min_digits`
  expect_equal(loo:::.se_digits(45), 1)
  # nothing usable is left
  expect_equal(loo:::.se_digits(c(0, NA, Inf)), 2)
  expect_equal(loo:::.se_digits(NULL), 2)
})

test_that(".measure_digits is fixed for a measure on a fixed scale", {
  expect_equal(loo:::.measure_digits("elpd"), 1)
  expect_equal(loo:::.measure_digits("ic"), 1)
  expect_equal(loo:::.measure_digits("mlpd"), 3)
  for (m in c("r2", "acc", "bacc", "brier")) {
    expect_equal(loo:::.measure_digits(m), 3)
  }
  # on the scale of the data the standard error decides
  expect_equal(loo:::.measure_digits("rmse", 0.00031), 4)
  expect_equal(loo:::.measure_digits("mae", 0.45), 2)
  # a custom measure has no entry, so it follows the standard error too
  expect_equal(loo:::.measure_digits("my_measure", 2.3), 1)
})

test_that(".resolve_digits honours the user's `digits`", {
  expect_equal(loo:::.resolve_digits(NULL, "r2", 0.006), 3)
  expect_equal(loo:::.resolve_digits(2, "r2", 0.006), 2)
  expect_equal(loo:::.resolve_digits(c(r2 = 5), "r2", 0.006), 5)
  # a measure the vector does not name keeps its default
  expect_equal(loo:::.resolve_digits(c(r2 = 5), "elpd", 1.4), 1)
})

test_that(".fr keeps a small value out of scientific notation", {
  expect_equal(loo:::.fr(0.00031, 4), "0.0003")
  expect_equal(loo:::.fr(c(0.00031, 0), 4), c("0.0003", "0.0000"))
})
