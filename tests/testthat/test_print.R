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
        measure = measure$measure,
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