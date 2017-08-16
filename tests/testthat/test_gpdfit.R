library(loo)

context("gpdfit")
test_that("gpdfit returns correct result", {
  set.seed(123)
  x <- rexp(100)
  gpdfit_val <- unlist(gpdfit(x))
  expect_equal_to_reference(gpdfit_val, "gpdfit.rds")
})
