# Extracted from test_compare.R:267

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "loo", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
set.seed(123)
LLarr <- example_loglik_array()
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- suppressWarnings(waic(LLarr))
w2 <- suppressWarnings(waic(LLarr2))

# test -------------------------------------------------------------------------
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
    rows <- out[grepl("^\\s+m[0-9]", out)]
    sub("^\\s*(\\S+).*$", "\\1", rows)
  }
for (comp in list(
    suppressMessages(model_compare(pms)),
    suppressMessages(model_compare(pms, rank_by = "mse"))
  )) {
    for (measure in c("elpd", "r2", "mse", "mae")) {
      diff_col <- if (measure == "elpd") "elpd_diff" else paste0(measure, "_diff")
      expected <- comp$model[order(comp[[diff_col]], decreasing = TRUE)]
      expect_equal(printed_order(comp, measure), expected)
      # the best model on the measure leads, and the table runs downhill
      expect_equal(expected[[1L]], comp$model[[which.max(comp[[diff_col]])]])
      expect_false(is.unsorted(rev(comp[[diff_col]][order(comp[[diff_col]])])))
    }
  }
