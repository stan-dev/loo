library(loo)

# test gpdfit -------------------------------------------------------------
context("gpdfit")
set.seed(123)
x <- rexp(100)
gpdfit_val <- unlist(gpdfit(x))
gpdfit_ans <- c(k = 0.02740303, sigma = 1.018298)
gpdfit_diff <- gpdfit_val - gpdfit_ans

test_that("gpdfit returns correct result", {
  expect_true(all(gpdfit_diff < 1e-6))
})
