library(loo)

# test gpdfit -------------------------------------------------------------
context("gpdfit")
set.seed(123)
x <- rexp(100)
gpdfit_val <- unlist(gpdfit(x))
gpdfit_ans <- structure(c(0.0274030348712631, 1.01829821712701),
                        .Names = c("k", "sigma"))

test_that("gpdfit returns correct result", {
  expect_equal(gpdfit_val, gpdfit_ans)
})
