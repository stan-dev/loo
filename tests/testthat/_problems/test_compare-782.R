# Extracted from test_compare.R:782

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

# test -------------------------------------------------------------------------
res <- readRDS("data-for-tests/test_data_roaches_compare.Rds")
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
    make(res$loo_p_m1, res$mupred_m1, res$ylp_m1, my_rmse),
    make(res$loo_p_m2, res$mupred_m2, res$ylp_m2, my_rmse)
  )
expect_equal(
    attr(pms[[1L]], "measure_compare_meta")$my_rmse$diff_method,
    "custom"
  )
expect_null(attr(pms[[1L]], "measure_compare_meta")$my_rmse$se_diff_fun)
expect_equal(
    loo:::.measure_pointwise_diff_method(pms, "my_rmse_loo"),
    "custom"
  )
expect_error(
    suppressMessages(model_compare(pms)),
    "my_rmse.*custom measure.*must be supplied"
  )
comp_null <- suppressMessages(model_compare(pms, custom_se_fn = NULL))
expect_false(is.na(comp_null$my_rmse_diff[[2L]]))
expect_true(all(is.na(comp_null$my_rmse_se_diff)))
expect_true(is.na(
    loo:::.pair_measure_stats(pms[[2L]], pms[[1L]], "my_rmse_loo", loos = pms)["se"]
  ))
comp_fn <- suppressMessages(
    model_compare(pms, custom_se_fn = loo:::.se_diff_rmse)
  )
expect_false(any(is.na(comp_fn$my_rmse_se_diff)))
expect_gt(comp_fn$my_rmse_se_diff[[2L]], 0)
