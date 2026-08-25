# Extracted from test_compare.R:900

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
expect_true(attr(declared$m1, "measure_compare_meta")$my_mse$loss)
expect_true(loo:::.measure_is_loss("my_mse_loo", declared))
expect_true(loo:::.measure_lower_is_better("my_mse_loo", declared))
expect_false(loo:::.measure_lower_is_better("my_mse_loo", plain))
expect_message(
    comp <- model_compare(declared, custom_se_fn = "mean"),
    "my_mse.*utility scale"
  )
comp_plain <- suppressMessages(model_compare(plain, custom_se_fn = "mean"))
expect_equal(attr(comp, "sign_converted_measures"), "my_mse")
expect_length(attr(comp_plain, "sign_converted_measures"), 0L)
expect_equal(comp$my_mse_diff, -comp_plain$my_mse_diff)
