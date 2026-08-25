# Extracted from test_compare.R:261

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
