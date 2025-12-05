library(loo)

context("generalized pareto")

test_that("gpdfit returns correct result", {
  set.seed(123)
  x <- rexp(100)
  gpdfit_val_old <- unlist(gpdfit(x, wip=FALSE, min_grid_pts = 80))
  expect_equal_to_reference(gpdfit_val_old, "reference-results/gpdfit_old.rds")

  gpdfit_val_wip <- unlist(gpdfit(x, wip=TRUE, min_grid_pts = 80))
  expect_equal_to_reference(gpdfit_val_wip, "reference-results/gpdfit.rds")

  gpdfit_val_wip_default_grid <- unlist(gpdfit(x, wip=TRUE))
  expect_equal_to_reference(gpdfit_val_wip_default_grid, "reference-results/gpdfit_default_grid.rds")
})

test_that("qgpd returns the correct result ", {
  probs <- seq(from = 0, to = 1, by = 0.25)
  q1 <- qgpd(probs, k = 1, sigma = 1)
  expect_equal(q1, c(0, 1/3, 1, 3, Inf))

  q2 <- qgpd(probs, k = 1, sigma = 0)
  expect_true(all(is.nan(q2)))
})
